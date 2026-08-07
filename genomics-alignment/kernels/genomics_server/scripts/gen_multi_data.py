#!/usr/bin/env python3
"""Generate multi-iteration data for mcast_array v3m.

Inputs  (NUM_ROWS files): input_r.txt  — one per AIE row input PLIO.
Goldens (NUM_COLS files): golden_c.txt — one per AIE column output PLIO.

golden_c.txt format per iteration: [count, score_row(R-1)_ref_c, ..., score_row0_ref_c, ...]
where ref indices stride by NUM_COLS: c, c+NUM_COLS, c+2*NUM_COLS, ...
This matches the column-reduce output order from TileUpperV3/TileBottomV3.

All kernels use 2-bit packing: 16 bases/int32, base j = (word >> j*2) & 0x03, A=0 C=1 G=2 T=3.
Kernels differ only in their scoring function:
  sw (default)  SW BLOSUM62 gap-affine local alignment (maximize score)
  wfa           gap-affine edit distance (minimize)
  bitap         Wu-Manber shift-OR edit distance, MAX_K=5 (minimize; -1 if no match)

Usage (synthetic):
  python gen_multi_data.py <out_dir> <num_rows> <num_cols> <refs_per_batch> <num_iters> <buf_size> [seed]
      [--kernel wfa|sw|bitap] [--read-len N] [--bitap-max-k N]
      [--snp-rate F] [--ins-rate F] [--del-rate F] [--error-rate F]

Usage (real Illumina reads):
  python gen_multi_data.py <out_dir> <num_rows> <num_cols> <refs_per_batch> <num_iters> <buf_size> [seed]
      [--kernel wfa|sw|bitap] --fastq reads.fastq [--ref-fasta ref.fa] [--read-len N] [--bitap-max-k N]
      [--snp-rate F] [--ins-rate F] [--del-rate F] [--error-rate F]

Usage (real BAM -- mapped reads at true mapping positions):
  python gen_multi_data.py <out_dir> <num_rows> <num_cols> <refs_per_batch> <num_iters> <buf_size> [seed]
      [--kernel wfa|sw|bitap] --bam reads.bam --ref-fasta ref.fa [--min-mapq N] [--read-len N] [--bitap-max-k N]
      [--snp-rate F] [--ins-rate F] [--del-rate F] [--error-rate F]

Mutation rates default to snp=0.01, ins=0.001, del=0.001, error=0.005 --
matching mutate_sequence() in scripts/generate_test_data.py (and
scripts/gen_interseq_data.py, which imports it), so synthetic data is
comparable across all three generators. In synthetic mode, near-position
decoys (shift +/-5bp) use 0.75x these rates; far-position decoys (+/-15bp)
use the full rates -- see mutate_sequence() call sites in _gen_row_synthetic.
"""
import sys, os, random, gzip, math, multiprocessing as mp, functools

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'scripts'))
from generate_test_data import wfa_align

try:
    import pysam
    _PYSAM_AVAILABLE = True
except ImportError:
    _PYSAM_AVAILABLE = False

BASE = 'ACGT'
BASE_TO_INT = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
BLOSUM62 = [[4, 0, 0, 0], [0, 9, -3, -1], [0, -3, 6, -2], [0, -1, -2, 5]]

# ── sequence utilities ────────────────────────────────────────────────────────

def _clean_bases(seq):
    """Replace non-ACGT characters (N, ambiguity codes) with a random valid base."""
    return ''.join(b if b in BASE_TO_INT else random.choice(BASE) for b in seq)

def random_dna(n):
    return ''.join(random.choice(BASE) for _ in range(n))

def mutate_sequence(seq, snp_rate=0.01, ins_rate=0.001, del_rate=0.001, error_rate=0.005):
    """Independent per-base SNP/insertion/deletion, plus a second-pass sequencing
    error rate. Matches mutate_sequence() in scripts/generate_test_data.py (and
    scripts/gen_interseq_data.py, which imports it directly) so synthetic data
    across all three generators uses the same mutation model."""
    result = []
    for base in seq:
        r = random.random()
        if r < del_rate:
            continue
        if r < del_rate + ins_rate:
            result.append(random.choice(BASE)); result.append(base); continue
        if r < del_rate + ins_rate + snp_rate:
            result.append(random.choice([b for b in BASE if b != base])); continue
        result.append(base)
    mutated = []
    for base in result:
        if random.random() < error_rate:
            mutated.append(random.choice([b for b in BASE if b != base]))
        else:
            mutated.append(base)
    return ''.join(mutated)

# ── packing (2-bit for both kernels) ─────────────────────────────────────────

def pack_bases_2bit(seq):
    """16 bases per int32, base at position j: (word >> (j*2)) & 0x03"""
    words = []
    for i in range(0, len(seq), 16):
        w = 0
        for j in range(min(16, len(seq) - i)):
            w |= BASE_TO_INT[seq[i + j]] << (j * 2)
        words.append(w)
    return words

# ── scoring ───────────────────────────────────────────────────────────────────

def bitap_edit_distance(pattern, text, max_k=5):
    """Wu-Manber shift-OR approximate matching. Matches bitap_kernel.cc exactly.
    Bit convention: 0=active (reachable), 1=inactive — same as C kernel.
    Returns min edit distance 0..max_k, or -1 if no alignment within max_k edits."""
    plen = len(pattern)
    if plen == 0:
        return 0
    k = min(max_k, plen)
    mask = (1 << plen) - 1

    # PM[base]: bit i=0 if pattern[i]==base (match), 1=mismatch
    PM = {b: mask for b in 'ACGT'}
    for i, b in enumerate(pattern):
        PM[b] &= ~(1 << i)

    # R[d]: first d bits cleared = d free deletions from start (matches kernel init)
    R = [mask] * (k + 1)
    for d in range(k + 1):
        for i in range(min(d, plen)):
            R[d] &= ~(1 << i)

    accept = 1 << (plen - 1)

    for ch in text:
        pm = PM.get(ch, mask)   # unknown base = all-mismatch
        old_R = R[:]

        # d=0: exact match — R[0] = (R[0]<<1) | pm
        R[0] = ((old_R[0] << 1) | pm) & mask

        # d=1..k: mirrors BITAP_LEVEL macro in bitap_kernel.cc:
        #   R[d] = ((old[d]<<1)|pm) & (old[d-1]<<1) & old[d-1] & (new[d-1]<<1)
        for d in range(1, k + 1):
            t  = ((old_R[d] << 1) | pm) & mask     # match extension from d
            t &= ((old_R[d-1] << 1) & mask)         # substitution  (old d-1, shifted)
            t &= old_R[d-1]                          # insertion     (old d-1)
            t &= ((R[d-1] << 1) & mask)             # deletion      (new d-1, shifted)
            R[d] = t

    # Final state: find min d where bit(plen-1) is clear (full pattern reached)
    for d in range(k + 1):
        if (R[d] & accept) == 0:
            return d
    return -1

def sw_blosum62(q, r):
    """Smith-Waterman local alignment, BLOSUM62 DNA, gap_open=5, gap_extend=1.
    E/F floored at 0 to match AIE saturating_subu behaviour."""
    GO, GE = 5, 1
    B = BASE_TO_INT
    m, n = len(q), len(r)
    H = [[0] * (n + 1) for _ in range(m + 1)]
    E = [[0] * (n + 1) for _ in range(m + 1)]
    F = [[0] * (n + 1) for _ in range(m + 1)]
    best = 0
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            sub = BLOSUM62[B[q[i-1]]][B[r[j-1]]]
            E[i][j] = max(H[i][j-1] - GO, E[i][j-1] - GE)
            F[i][j] = max(H[i-1][j] - GO, F[i-1][j] - GE)
            H[i][j] = max(H[i-1][j-1] + sub, E[i][j], F[i][j], 0)
            best = max(best, H[i][j])
    return best

# ── real Illumina data readers ────────────────────────────────────────────────

def _open_maybe_gzip(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')

def read_fasta(path):
    """Parse FASTA (plain or .gz). Returns {contig_name: sequence_str}."""
    sequences = {}
    name, parts = None, []
    with _open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    sequences[name] = _clean_bases(''.join(parts))
                name  = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.upper())
    if name is not None:
        sequences[name] = _clean_bases(''.join(parts))
    print(f'  FASTA: loaded {len(sequences)} contig(s) from {path}')
    return sequences

def read_fastq(path, max_reads=0, min_len=0, max_len=0):
    """Parse FASTQ (plain or .gz). Returns list of (read_name, sequence)."""
    reads = []
    with _open_maybe_gzip(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip().upper()
            fh.readline()   # +
            fh.readline()   # quality
            if not header.startswith('@'):
                continue
            seq = _clean_bases(seq)
            if min_len and len(seq) < min_len:
                continue
            if max_len and len(seq) > max_len:
                seq = seq[:max_len]
            reads.append((header[1:].split()[0], seq))
            if max_reads and len(reads) >= max_reads:
                break
    print(f'  FASTQ: loaded {len(reads)} reads from {path}')
    return reads

def read_bam_pairs(bam_path, fasta_seqs, max_pairs, read_len, min_mapq=20):
    """Extract (read, ref_window) pairs from a mapped BAM at true mapping positions."""
    if not _PYSAM_AVAILABLE:
        print('ERROR: pysam is required for --bam mode. Install with: pip install pysam')
        sys.exit(1)
    pairs = []
    bam = pysam.AlignmentFile(bam_path, 'rb')
    skipped = 0
    for read in bam.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            skipped += 1
            continue
        if read.mapping_quality < min_mapq:
            skipped += 1
            continue
        chrom = read.reference_name
        pos   = read.reference_start          # 0-based
        if chrom not in fasta_seqs:
            skipped += 1
            continue
        ref = fasta_seqs[chrom]
        if pos + read_len > len(ref):
            skipped += 1
            continue
        raw_seq = read.query_sequence
        if raw_seq is None or len(raw_seq) < read_len:
            skipped += 1
            continue
        seq   = _clean_bases(raw_seq[:read_len])   # real Illumina read
        ref_window = ref[pos:pos + read_len]        # reference window at true mapping site
        pairs.append((seq, ref_window))
        if len(pairs) >= max_pairs:
            break
    bam.close()
    print(f'  BAM: {len(pairs)} pairs selected, {skipped} reads skipped '
          f'(unmapped/secondary/low-mapq/short)')
    if len(pairs) < max_pairs:
        print(f'  Warning: only {len(pairs)} usable reads found (requested {max_pairs})')
    return pairs

def fasta_window(sequences, length):
    """Random window of `length` bases from FASTA, weighted by contig length."""
    names   = list(sequences.keys())
    lengths = [len(sequences[n]) for n in names]
    r       = random.random() * sum(lengths)
    cumsum  = 0
    chosen  = names[0]
    for n, l in zip(names, lengths):
        cumsum += l
        if r < cumsum:
            chosen = n
            break
    ref_seq = sequences[chosen]
    max_pos = len(ref_seq) - length
    if max_pos <= 0:
        return ref_seq[:length]
    pos = random.randint(0, max_pos)
    return ref_seq[pos:pos + length]

# ── shared helpers ────────────────────────────────────────────────────────────

def build_batch_words(query, refs, buf_size):
    """Encode one batch into exactly buf_size int32 words (2-bit packing)."""
    words = [len(refs), len(query)] + pack_bases_2bit(query)
    for ref in refs:
        words.append(len(ref))
        words.extend(pack_bases_2bit(ref))
    while len(words) < buf_size:
        words.append(0)
    return words[:buf_size]

def write_input(path, all_words):
    with open(path, 'w') as f:
        for w in all_words:
            if w >= (1 << 31): w -= (1 << 32)
            f.write(f'{w}\n')

# ── per-column golden writer (top-level for Pool pickling) ────────────────────

def _write_col_golden(args):
    c, out_dir, col_scores, num_rows, num_iters = args
    # col_scores[r][it] = list of scores for this column's refs
    with open(os.path.join(out_dir, f'golden_{c}.txt'), 'w') as f:
        for it in range(num_iters):
            count = num_rows * len(col_scores[0][it])
            f.write(f'{count}\n')
            for r in range(num_rows - 1, -1, -1):   # bottom row first
                for s in col_scores[r][it]:
                    f.write(f'{s}\n')

def write_column_goldens(out_dir, all_scores, num_rows, num_cols, num_iters, refs_per_batch):
    """Write golden_c.txt for each column c in parallel."""
    col_refs_cache = [list(range(c, refs_per_batch, num_cols)) for c in range(num_cols)]
    # Pre-extract per-column scores so each worker receives only what it needs.
    tasks = []
    for c in range(num_cols):
        col_refs = col_refs_cache[c]
        col_scores = [[[all_scores[r][it][i] for i in col_refs]
                       for it in range(num_iters)]
                      for r in range(num_rows)]
        tasks.append((c, out_dir, col_scores, num_rows, num_iters))
    nw = min(num_cols, mp.cpu_count())
    with mp.Pool(nw) as pool:
        pool.map(_write_col_golden, tasks)

# ── row workers (top-level for Pool pickling) ─────────────────────────────────

def _gen_row_synthetic(args):
    (r, num_iters, refs_per_batch, buf_size, QLEN, RLEN,
     genome, kernel, seed, max_k, rates) = args
    snp_rate, ins_rate, del_rate, error_rate = rates
    random.seed(seed + r * 1337)
    score = (sw_blosum62                                           if kernel == 'sw' else
             wfa_align                                             if kernel == 'wfa'         else
             functools.partial(bitap_edit_distance, max_k=max_k))
    GENOME_LEN = len(genome)
    all_words, row_scores = [], []
    for it in range(num_iters):
        true_pos = random.randint(0, GENOME_LEN - QLEN)
        query    = genome[true_pos:true_pos + QLEN]
        refs, scores = [], []
        for ri in range(refs_per_batch):
            if ri == 0:
                ref = mutate_sequence(query, snp_rate, ins_rate, del_rate, error_rate)
            elif ri < refs_per_batch * 3 // 4:
                # near-position decoy: 0.75x the base rates (closer position = more similar,
                # harder to distinguish from the true match -- see conversation for why this
                # tiering exists: avoids the degenerate all-decoys-equally-wrong distribution).
                shift = random.randint(-5, 5)
                pos   = max(0, min(GENOME_LEN - RLEN, true_pos + shift))
                ref   = mutate_sequence(genome[pos:pos + RLEN],
                                         snp_rate * 0.75, ins_rate * 0.75, del_rate * 0.75, error_rate)
            else:
                shift = random.randint(-15, 15)
                pos   = max(0, min(GENOME_LEN - RLEN, true_pos + shift))
                ref   = mutate_sequence(genome[pos:pos + RLEN], snp_rate, ins_rate, del_rate, error_rate)
            refs.append(ref)
            scores.append(score(query, ref))
        row_scores.append(scores)
        all_words.extend(build_batch_words(query, refs, buf_size))
    return r, all_words, row_scores

def _gen_row_illumina(args):
    (r, num_iters, refs_per_batch, buf_size, QLEN, RLEN,
     fastq_reads, fasta_seqs, kernel, seed, read_start, max_k, rates) = args
    snp_rate, ins_rate, del_rate, error_rate = rates
    random.seed(seed + r * 1337)
    score = (sw_blosum62                                           if kernel == 'sw' else
             wfa_align                                             if kernel == 'wfa'         else
             functools.partial(bitap_edit_distance, max_k=max_k))
    all_words, row_scores = [], []
    for it in range(num_iters):
        _, seq = fastq_reads[(read_start + it) % len(fastq_reads)]
        query = seq[:QLEN] if len(seq) >= QLEN else seq + random_dna(QLEN - len(seq))
        refs, scores = [], []
        for _ in range(refs_per_batch):
            ref = (fasta_window(fasta_seqs, RLEN) if fasta_seqs
                   else mutate_sequence(query, snp_rate, ins_rate, del_rate, error_rate))
            refs.append(ref)
            scores.append(score(query, ref))
        row_scores.append(scores)
        all_words.extend(build_batch_words(query, refs, buf_size))
    return r, all_words, row_scores

def _gen_row_bam(args):
    (r, num_iters, refs_per_batch, buf_size, QLEN, RLEN,
     bam_pairs, kernel, seed, pair_start, max_k, rates) = args
    snp_rate, ins_rate, del_rate, error_rate = rates
    random.seed(seed + r * 1337)
    score = (sw_blosum62                                           if kernel == 'sw' else
             wfa_align                                             if kernel == 'wfa'         else
             functools.partial(bitap_edit_distance, max_k=max_k))
    all_words, row_scores = [], []
    for it in range(num_iters):
        seq, true_ref = bam_pairs[(pair_start + it) % len(bam_pairs)]
        query = seq[:QLEN] if len(seq) >= QLEN else seq + random_dna(QLEN - len(seq))
        # true mapping-site reference window (ground truth)
        true_window = true_ref[:RLEN] if len(true_ref) >= RLEN else true_ref + random_dna(RLEN - len(true_ref))
        refs, scores = [], []
        for ri in range(refs_per_batch):
            # ri==0 is the exact true window; the rest are near-variants of it
            # (not unrelated genome windows) so scores stay within max_k instead of -1.
            ref = (true_window if ri == 0
                   else mutate_sequence(true_window, snp_rate, ins_rate, del_rate, error_rate))
            refs.append(ref)
            scores.append(score(query, ref))
        row_scores.append(scores)
        all_words.extend(build_batch_words(query, refs, buf_size))
    return r, all_words, row_scores

# ── argument parsing ──────────────────────────────────────────────────────────

def parse_args():
    positional, flags = [], {}
    i = 1
    while i < len(sys.argv):
        a = sys.argv[i]
        if a.startswith('--'):
            key = a.lstrip('-')
            if '=' in key:
                k, v = key.split('=', 1)
                flags[k] = v
            elif i + 1 < len(sys.argv) and not sys.argv[i + 1].startswith('--'):
                flags[key] = sys.argv[i + 1]
                i += 1
            else:
                flags[key] = True
        else:
            positional.append(a)
        i += 1
    return positional, flags

# ── main ──────────────────────────────────────────────────────────────────────

def main():
    positional, flags = parse_args()

    if len(positional) < 6:
        print(__doc__)
        sys.exit(1)

    out_dir        = positional[0]
    num_rows       = int(positional[1])
    num_cols       = int(positional[2])
    refs_per_batch = int(positional[3])
    num_iters      = int(positional[4])
    buf_size       = int(positional[5])
    seed           = int(positional[6]) if len(positional) > 6 else 42

    kernel    = flags.get('kernel', 'sw')
    read_len  = int(flags.get('read-len', 128))
    max_k     = int(flags.get('bitap-max-k', 5))
    fastq     = flags.get('fastq')
    bam       = flags.get('bam')
    ref_fasta = flags.get('ref-fasta')
    min_mapq  = int(flags.get('min-mapq', 20))

    # Matches generate_test_data.py's mutate_sequence() defaults, so synthetic
    # data across all three generators (this file, generate_test_data.py,
    # gen_interseq_data.py) uses the same mutation model.
    snp_rate   = float(flags.get('snp-rate', 0.02))
    ins_rate   = float(flags.get('ins-rate', 0.0001))
    del_rate   = float(flags.get('del-rate', 0.0001))
    error_rate = float(flags.get('error-rate', 0.005))
    rates = (snp_rate, ins_rate, del_rate, error_rate)

    if kernel not in ('sw', 'wfa', 'bitap'):
        print(f'ERROR: unknown kernel "{kernel}". Choose: sw | wfa | bitap')
        sys.exit(1)

    if bam and not ref_fasta:
        print('ERROR: --bam requires --ref-fasta')
        sys.exit(1)

    if read_len > 128:
        print(f'WARNING: --read-len {read_len} > MAX_SEQ_LEN (128), clamping to 128.')
        read_len = 128

    random.seed(seed)
    os.makedirs(out_dir, exist_ok=True)

    QLEN = RLEN = read_len
    nworkers = min(num_rows, mp.cpu_count())

    # all_scores[r][it] = list of refs_per_batch scores
    all_scores = [None] * num_rows

    # ── real BAM (mapped reads at true mapping positions) ────────────────────
    if bam:
        print(f'Mode: REAL BAM -- {bam}')
        fasta_seqs = read_fasta(ref_fasta)
        bam_pairs = read_bam_pairs(bam, fasta_seqs,
                                   max_pairs=num_rows * num_iters * 2,
                                   read_len=QLEN, min_mapq=min_mapq)
        if not bam_pairs:
            print('ERROR: no usable pairs from BAM'); sys.exit(1)

        total_pairs = num_rows * num_iters * refs_per_batch
        est_sec = total_pairs * (QLEN * RLEN) / 5e6
        print(f'Kernel: {kernel}  packing: 2-bit (16/word)  read_len: {QLEN}bp  '
              f'refs/batch: {refs_per_batch}  num_cols: {num_cols}')
        print(f'Scoring {total_pairs} pairs ({QLEN}x{RLEN} DP) using {nworkers} workers'
              f' -- est. {math.ceil(est_sec / nworkers)}s ...')

        tasks = [
            (r, num_iters, refs_per_batch, buf_size, QLEN, RLEN,
             bam_pairs, kernel, seed, r * num_iters, max_k, rates)
            for r in range(num_rows)
        ]
        with mp.Pool(nworkers) as pool:
            for r, all_words, row_scores in pool.imap_unordered(_gen_row_bam, tasks):
                all_scores[r] = row_scores
                write_input(os.path.join(out_dir, f'input_{r}.txt'), all_words)
                print(f'  row {r}: {num_iters} iters x {refs_per_batch} refs  done')

    # ── real Illumina reads ───────────────────────────────────────────────────
    elif fastq:
        print(f'Mode: REAL Illumina -- {fastq}')
        fastq_reads = read_fastq(fastq,
                                 max_reads=num_rows * num_iters * 2,
                                 min_len=20,
                                 max_len=QLEN)
        if not fastq_reads:
            print('ERROR: no reads loaded from FASTQ'); sys.exit(1)
        fasta_seqs = read_fasta(ref_fasta) if ref_fasta else None

        total_pairs = num_rows * num_iters * refs_per_batch
        est_sec = total_pairs * (QLEN * RLEN) / 5e6
        print(f'Kernel: {kernel}  packing: 2-bit (16/word)  read_len: {QLEN}bp  '
              f'refs/batch: {refs_per_batch}  num_cols: {num_cols}')
        print(f'Scoring {total_pairs} pairs ({QLEN}x{RLEN} DP) using {nworkers} workers'
              f' -- est. {math.ceil(est_sec / nworkers)}s ...')

        tasks = [
            (r, num_iters, refs_per_batch, buf_size, QLEN, RLEN,
             fastq_reads, fasta_seqs, kernel, seed, r * num_iters, max_k, rates)
            for r in range(num_rows)
        ]
        with mp.Pool(nworkers) as pool:
            for r, all_words, row_scores in pool.imap_unordered(_gen_row_illumina, tasks):
                all_scores[r] = row_scores
                write_input(os.path.join(out_dir, f'input_{r}.txt'), all_words)
                print(f'  row {r}: {num_iters} iters x {refs_per_batch} refs  done')

    # ── synthetic ─────────────────────────────────────────────────────────────
    else:
        GENOME_LEN = 10000
        genome = random_dna(GENOME_LEN)
        total_pairs = num_rows * num_iters * refs_per_batch
        est_sec = total_pairs * (QLEN * RLEN) / 5e6
        print(f'Mode: synthetic  kernel: {kernel}  packing: 2-bit (16/word)  '
              f'read_len: {QLEN}bp  num_cols: {num_cols}')
        print(f'Scoring {total_pairs} pairs ({QLEN}x{RLEN} DP) using {nworkers} workers'
              f' -- est. {math.ceil(est_sec / nworkers)}s ...')

        tasks = [
            (r, num_iters, refs_per_batch, buf_size, QLEN, RLEN,
             genome, kernel, seed, max_k, rates)
            for r in range(num_rows)
        ]
        with mp.Pool(nworkers) as pool:
            for r, all_words, row_scores in pool.imap_unordered(_gen_row_synthetic, tasks):
                all_scores[r] = row_scores
                write_input(os.path.join(out_dir, f'input_{r}.txt'), all_words)
                print(f'  row {r}: {num_iters} iters x {refs_per_batch} refs  done')

    # ── write per-column goldens in parallel (both modes) ─────────────────────
    print(f'Writing {num_cols} golden files using {min(num_cols, mp.cpu_count())} workers ...')
    write_column_goldens(out_dir, all_scores, num_rows, num_cols, num_iters, refs_per_batch)

    total = num_rows * num_iters * refs_per_batch
    print(f'\nWrote {num_rows} input files, {num_cols} golden files')
    print(f'Total: {num_rows} rows x {num_iters} iters x {refs_per_batch} refs = {total} refs')

if __name__ == '__main__':
    main()
