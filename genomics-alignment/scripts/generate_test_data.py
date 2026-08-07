#!/usr/bin/env python3
"""
Realistic test data generator for short read alignment kernels.

Generates test datasets that mirror real Illumina sequencing data:
- Reads of configurable length (75-250 bp, default 100/150)
- Controlled mutation profiles: SNPs, insertions, deletions, sequencing errors
- Reference sequences from synthetic DNA or real chromosome fragments
- Ground-truth alignment scores for all algorithm variants
- Outputs in both human-readable and AIE PLIO-compatible formats

Scoring: BLOSUM62 DNA substitution matrix (4x4, A/C/G/T)
  Match +2, Transition (A<->G, C<->T) -1, Transversion -2
  Must match BLOSUM62_DNA[][] in kernels/common/dna_utils.hh

Synthetic mode (default):
    python generate_test_data.py --preset illumina --out-dir data/

Real Illumina mode (requires FASTQ + FASTA):
    python generate_test_data.py --fastq reads.fastq --ref-fasta ref.fa \
        --num-pairs 1000 --out-dir data/illumina_real/

    # PhiX174 (Illumina control) example:
    # 1. Download reference: https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1?report=fasta
    # 2. Download reads:     https://www.ebi.ac.uk/ena/browser/view/SRR000001
    python generate_test_data.py --fastq phix_reads.fastq \
        --ref-fasta phix174.fa --num-pairs 500 --out-dir data/phix/
"""

import argparse
import gzip
import json
import os
import random
import sys
from dataclasses import dataclass
from typing import List, Tuple, Optional

try:
    import pysam
    _PYSAM_AVAILABLE = True
except ImportError:
    _PYSAM_AVAILABLE = False

BASES = "ACGT"
BASE_TO_INT = {"A": 0, "C": 1, "G": 2, "T": 3}
INT_TO_BASE = {0: "A", 1: "C", 2: "G", 3: "T"}

BLOSUM62_DNA = [
#   A   C   G   T
    [ 4, 0, 0, 0],   # A: match=4, transv C=0, transit G=0, transv T=0
    [0,  9, -3, -1],   # C: transv A=0, match=9, transv G=-3, transit T=-1
    [0, -3,  6, -2],   # G: transit A=0, transv C=-3, match=6, transv T=-2
    [0, -1, -2,  5],    # T: transv A=0, transit C=-1, transv G=-2, match=5
]

def _blosum62(b1: str, b2: str) -> int:
    """Lookup BLOSUM62 score for two DNA bases."""
    return BLOSUM62_DNA[BASE_TO_INT[b1]][BASE_TO_INT[b2]]

# ============================================================================
# Reference Implementations (ground truth)
# ============================================================================

def smith_waterman(seq1: str, seq2: str,
                   gap_open=-5, gap_extend=-1) -> int:
    """Local alignment — BLOSUM62 DNA substitution, affine gaps (corrected Gotoh)."""
    m, n = len(seq1), len(seq2)
    M = [[0] * (n + 1) for _ in range(m + 1)]
    I = [[0] * (n + 1) for _ in range(m + 1)]
    D = [[0] * (n + 1) for _ in range(m + 1)]
    best = 0
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = _blosum62(seq1[i-1], seq2[j-1])
            I[i][j] = max(M[i-1][j] + gap_open,  I[i-1][j] + gap_extend)
            D[i][j] = max(M[i][j-1] + gap_open,  D[i][j-1] + gap_extend)
            M[i][j] = max(0, M[i-1][j-1] + s, I[i][j], D[i][j])
            best = max(best, M[i][j])
    return best


def needleman_wunsch(seq1: str, seq2: str,
                     gap_open=-5, gap_extend=-1) -> int:
    """Global alignment — BLOSUM62 DNA substitution, affine gaps (corrected Gotoh)."""
    m, n = len(seq1), len(seq2)
    NEG_INF = -999999
    M = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    I = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    D = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    M[0][0] = 0
    for i in range(1, m + 1):
        I[i][0] = gap_open + (i - 1) * gap_extend
        M[i][0] = I[i][0]
    for j in range(1, n + 1):
        D[0][j] = gap_open + (j - 1) * gap_extend
        M[0][j] = D[0][j]
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = _blosum62(seq1[i-1], seq2[j-1])
            I[i][j] = max(M[i-1][j] + gap_open,  I[i-1][j] + gap_extend)
            D[i][j] = max(M[i][j-1] + gap_open,  D[i][j-1] + gap_extend)
            M[i][j] = max(M[i-1][j-1] + s, I[i][j], D[i][j])
    return max(M[m][n], I[m][n], D[m][n])


def banded_smith_waterman(seq1: str, seq2: str, bandwidth=16,
                          gap_open=-5, gap_extend=-1) -> int:
    """Banded local alignment — BLOSUM62 DNA substitution, corrected Gotoh."""
    m, n = len(seq1), len(seq2)
    M = [[0] * (n + 1) for _ in range(m + 1)]
    I = [[0] * (n + 1) for _ in range(m + 1)]
    D = [[0] * (n + 1) for _ in range(m + 1)]
    best = 0
    for i in range(1, m + 1):
        j_start = max(1, i - bandwidth)
        j_end   = min(n, i + bandwidth)
        for j in range(j_start, j_end + 1):
            s = _blosum62(seq1[i-1], seq2[j-1])
            I[i][j] = max(M[i-1][j] + gap_open,  I[i-1][j] + gap_extend)
            D[i][j] = max(M[i][j-1] + gap_open,  D[i][j-1] + gap_extend)
            M[i][j] = max(0, M[i-1][j-1] + s, I[i][j], D[i][j])
            best = max(best, M[i][j])
    return best


def bitap_edit_distance(pattern: str, text: str, max_k: int = 5) -> int:
    """Wu-Manber shift-OR approximate matching. Matches bitap_kernel.cc exactly.
    Supports patterns up to 128 bp (BV_WORDS=4, 4×32-bit bitvectors).
    Bit convention: 0=active (reachable), 1=inactive.
    Returns min edit distance 0..max_k, or -1 if no alignment found within max_k."""
    plen = len(pattern)
    if plen == 0:
        return 0
    k = min(max_k, plen)
    mask = (1 << plen) - 1

    # PM[base]: bit i=0 if pattern[i]==base (match), 1=mismatch
    PM = {b: mask for b in BASES}
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


def myers_edit_distance(pattern: str, text: str) -> int:
    """Myers bit-vector edit distance (pattern <= 32 bases)."""
    m = len(pattern)
    assert m <= 32
    PM = [0] * 4
    for i, base in enumerate(pattern):
        PM[BASE_TO_INT[base]] |= (1 << i)
    VP, VN, score, mask = (1 << m) - 1, 0, m, 1 << (m - 1)
    for j in range(len(text)):
        Eq = PM[BASE_TO_INT[text[j]]]
        X  = Eq | VN
        D0 = ((VP + (X & VP)) ^ VP) | X
        HP = VN | ~(D0 | VP)
        HN = D0 & VP
        if HP & mask: score += 1
        if HN & mask: score -= 1
        HP = ((HP << 1) | 1) & ((1 << m) - 1 | (1 << m))
        HN = (HN << 1)
        VP = (HN | ~(D0 | HP)) & ((1 << m) - 1)
        VN = (D0 & HP)         & ((1 << m) - 1)
    return score


def edit_distance_dp(seq1: str, seq2: str) -> int:
    m, n = len(seq1), len(seq2)
    prev = list(range(n + 1))
    curr = [0] * (n + 1)
    for i in range(1, m + 1):
        curr[0] = i
        for j in range(1, n + 1):
            cost = 0 if seq1[i-1] == seq2[j-1] else 1
            curr[j] = min(prev[j-1] + cost, prev[j] + 1, curr[j-1] + 1)
        prev, curr = curr, prev
    return prev[n]


def bitparallel_lcs(seq1: str, seq2: str) -> int:
    m, n = len(seq1), len(seq2)
    prev = [0] * (n + 1)
    for i in range(1, m + 1):
        curr = [0] * (n + 1)
        for j in range(1, n + 1):
            curr[j] = prev[j-1] + 1 if seq1[i-1] == seq2[j-1] else max(prev[j], curr[j-1])
        prev = curr
    return prev[n]


def wfa_align(seq1: str, seq2: str, gap_open=5, gap_extend=1,
              mismatch_penalty=4, max_score=60) -> int:
    """Wavefront alignment with fixed mismatch penalty.

    Any mismatched base pair costs mismatch_penalty (default 4).
    gap_open and gap_extend are WFA penalties (positive, cost to open/extend a gap).
    A gap of length L costs gap_open + gap_extend*L (matches WFA2-lib's
    gap_open1 = score - gap_opening1 - gap_extension1), not gap_open + gap_extend*(L-1).
    Returns minimum alignment penalty (0 = identical sequences), or -1 if it
    exceeds max_score -- mirrors wfa_kernel.cc's fixed MAX_SCORE=60 search-depth
    cap (the hardware gives up and returns -1 beyond that, so golden data must
    use the same cutoff rather than searching to a much larger safe bound).
    """
    m, n = len(seq1), len(seq2)

    target_k      = m - n
    target_offset = m
    # Cap the search at the caller's max_score (default 60, matching the
    # kernel's hardware search-depth limit) -- but never search further than
    # the mathematically safe upper bound for tiny sequences.
    max_score     = min(max_score, (m + n) * mismatch_penalty)

    wf_m = [None] * (max_score + 1)
    wf_i = [None] * (max_score + 1)
    wf_d = [None] * (max_score + 1)
    wf_m[0] = {0: 0}

    def extend(wf_s):
        if wf_s is None:
            return
        for k in list(wf_s.keys()):
            off = wf_s[k]
            i, j = off, off - k
            while 0 <= i < m and 0 <= j < n and seq1[i] == seq2[j]:
                off += 1; i += 1; j += 1
            wf_s[k] = off

    def gv(wf, k):
        return -999999 if (wf is None) else wf.get(k, -999999)

    extend(wf_m[0])
    if gv(wf_m[0], target_k) >= target_offset:
        return 0

    for s in range(1, max_score + 1):
        s_open = s - gap_open - gap_extend

        # Insertions (gap in seq1)
        wf_i[s] = {}
        if s_open >= 0 and wf_m[s_open] is not None:
            for k_src, off in wf_m[s_open].items():
                k = k_src - 1
                wf_i[s][k] = max(wf_i[s].get(k, -999999), off + 1)
        if s >= gap_extend and wf_i[s - gap_extend]:
            for k_src, off in wf_i[s - gap_extend].items():
                k = k_src - 1
                wf_i[s][k] = max(wf_i[s].get(k, -999999), off + 1)

        # Deletions (gap in seq2)
        wf_d[s] = {}
        if s_open >= 0 and wf_m[s_open] is not None:
            for k_src, off in wf_m[s_open].items():
                k = k_src + 1
                wf_d[s][k] = max(wf_d[s].get(k, -999999), off)
        if s >= gap_extend and wf_d[s - gap_extend]:
            for k_src, off in wf_d[s - gap_extend].items():
                k = k_src + 1
                wf_d[s][k] = max(wf_d[s].get(k, -999999), off)

        # Substitutions — fixed mismatch penalty
        wf_m[s] = {}
        if s >= mismatch_penalty and wf_m[s - mismatch_penalty] is not None:
            for k, off in wf_m[s - mismatch_penalty].items():
                i, j = off, off - k
                if 0 <= i < m and 0 <= j < n and seq1[i] != seq2[j]:
                    wf_m[s][k] = max(wf_m[s].get(k, -999999), off + 1)

        # Merge insertions and deletions into match wavefront
        for k in (wf_i[s] or {}):
            wf_m[s][k] = max(wf_m[s].get(k, -999999), wf_i[s][k])
        for k in (wf_d[s] or {}):
            wf_m[s][k] = max(wf_m[s].get(k, -999999), wf_d[s][k])

        extend(wf_m[s])
        if gv(wf_m[s], target_k) >= target_offset:
            return s

    return -1



# ============================================================================
# Real Illumina data readers
# ============================================================================

def _open_maybe_gzip(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_fasta(path: str) -> dict:
    """Parse a FASTA file. Returns {seq_name: sequence_str}.
    Non-ACGT characters replaced with a random ACGT base.
    """
    sequences = {}
    name, parts = None, []
    with _open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = _clean_bases(''.join(parts))
                name  = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.upper())
    if name is not None:
        sequences[name] = _clean_bases(''.join(parts))
    print(f"  FASTA: loaded {len(sequences)} contig(s) from {path}")
    return sequences


def read_fastq(path: str, max_reads: int = 0,
               min_len: int = 0, max_len: int = 0) -> List[Tuple[str, str]]:
    """Parse a FASTQ file (plain or .gz). Returns list of (read_name, sequence)."""
    reads = []
    with _open_maybe_gzip(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq  = fh.readline().rstrip().upper()
            fh.readline()
            fh.readline()
            if not header.startswith("@"):
                continue
            seq = _clean_bases(seq)
            length = len(seq)
            if min_len and length < min_len:   continue
            if max_len and length > max_len:   seq = seq[:max_len]
            reads.append((header[1:].split()[0], seq))
            if max_reads and len(reads) >= max_reads:
                break
    print(f"  FASTQ: loaded {len(reads)} reads from {path}")
    return reads


def _clean_bases(seq: str) -> str:
    result = []
    for b in seq:
        result.append(b if b in BASE_TO_INT else random.choice(BASES))
    return ''.join(result)


def fasta_window(sequences: dict, length: int) -> Tuple[str, str, int]:
    """Extract a random window of `length` bases, weighted by contig length."""
    names   = list(sequences.keys())
    lengths = [len(sequences[n]) for n in names]
    total   = sum(lengths)
    r       = random.random() * total
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
        return chosen, ref_seq[:length], 0
    pos = random.randint(0, max_pos)
    return chosen, ref_seq[pos:pos + length], pos


# ============================================================================
# Synthetic sequence generation
# ============================================================================

def random_dna(length: int) -> str:
    return ''.join(random.choice(BASES) for _ in range(length))


def mutate_sequence(seq: str, snp_rate=0.01, ins_rate=0.001, del_rate=0.001,
                    error_rate=0.005) -> Tuple[str, dict]:
    result = []
    mutations = {"snps": 0, "insertions": 0, "deletions": 0, "errors": 0}
    for base in seq:
        r = random.random()
        if r < del_rate:
            mutations["deletions"] += 1; continue
        if r < del_rate + ins_rate:
            mutations["insertions"] += 1
            result.append(random.choice(BASES)); result.append(base); continue
        if r < del_rate + ins_rate + snp_rate:
            mutations["snps"] += 1
            result.append(random.choice([b for b in BASES if b != base])); continue
        result.append(base)
    mutated = []
    for base in result:
        if random.random() < error_rate:
            mutations["errors"] += 1
            mutated.append(random.choice([b for b in BASES if b != base]))
        else:
            mutated.append(base)
    return ''.join(mutated), mutations


def extract_read(reference: str, read_len: int) -> Tuple[str, int]:
    max_pos = len(reference) - read_len
    if max_pos <= 0: return reference[:read_len], 0
    pos = random.randint(0, max_pos)
    return reference[pos:pos + read_len], pos


# ============================================================================
# PLIO Format Packing (2-bit encoding: 16 bases per int32)
# ============================================================================

def pack_bases_to_int32(bases: str) -> List[int]:
    """Pack bases into int32 words (16 bases per word, 2 bits per base).
    Encoding A=0,C=1,G=2,T=3. Matches unpack_bases_2bit() in swg_kernel.cc.
    """
    words = []
    for i in range(0, len(bases), 16):
        word = 0
        for j in range(16):
            if i + j < len(bases):
                word |= BASE_TO_INT[bases[i + j]] << (j * 2)
        words.append(word)
    return words


def format_plio_pair(seq1: str, seq2: str) -> Tuple[List[int], List[int]]:
    """Format a pair for AIE PLIO input (2-bit encoding, 16 bases/int32)."""
    lens_word = (len(seq1) << 16) | len(seq2)
    return ([lens_word] + pack_bases_to_int32(seq1),
            [lens_word] + pack_bases_to_int32(seq2))


# ============================================================================
# Dataset Generation
# ============================================================================

@dataclass
class TestPreset:
    name: str
    ref_len: int
    read_len: int
    num_pairs: int
    snp_rate: float
    ins_rate: float
    del_rate: float
    error_rate: float
    description: str

PRESETS = {
    "tiny":     TestPreset("tiny",     64,    16,   10,  0.0,  0.0,   0.0,   0.0,
        "Identical sequences, 8-16bp. Unit tests."),
    "small":    TestPreset("small",    200,   32,   50,  0.02, 0.0,   0.0,   0.005,
        "Small with SNPs only. x86sim."),
    "medium":   TestPreset("medium",   1000,  100,  200, 0.01, 0.001, 0.001, 0.005,
        "100bp synthetic reads. aiesim."),
    "illumina": TestPreset("illumina", 5000,  150,  1000,0.01, 0.002, 0.002, 0.005,
        "Synthetic 150bp. Use --fastq + --ref-fasta for real Illumina data."),
    "stress":   TestPreset("stress",   10000, 250,  5000,0.02, 0.003, 0.003, 0.01,
        "Long reads, high mutation. Performance testing."),
}


def _build_pairs_synthetic(preset, reference):
    pairs = []
    for i in range(preset.num_pairs):
        read, pos = extract_read(reference, preset.read_len)
        mutated, mutations = mutate_sequence(
            read, preset.snp_rate, preset.ins_rate, preset.del_rate, preset.error_rate)
        pairs.append({
            "index": i, "source": "synthetic",
            "ref_pos": pos,
            "seq1": mutated,
            "seq2": read,
            "mutations": mutations,
        })
    return pairs


def _build_pairs_real(fastq_reads, fasta_seqs, num_pairs, read_len):
    pairs = []
    n = min(num_pairs, len(fastq_reads))
    if n < num_pairs:
        print(f"  Warning: only {n} reads available (requested {num_pairs})")
    indices = list(range(len(fastq_reads)))
    random.shuffle(indices)
    for i, idx in enumerate(indices[:n]):
        name, read_seq = fastq_reads[idx]
        seq1 = read_seq[:read_len]
        _, seq2, ref_pos = fasta_window(fasta_seqs, len(seq1))
        pairs.append({
            "index": i, "source": "real_illumina",
            "read_name": name,
            "ref_pos": ref_pos,
            "seq1": seq1,
            "seq2": seq2,
            "mutations": {},
        })
    return pairs


def _build_pairs_fastq_only(fastq_reads, num_pairs, read_len,
                             snp_rate, ins_rate, del_rate, error_rate):
    pairs = []
    n = min(num_pairs, len(fastq_reads))
    if n < num_pairs:
        print(f"  Warning: only {n} reads available (requested {num_pairs})")
    indices = list(range(len(fastq_reads)))
    random.shuffle(indices)
    for i, idx in enumerate(indices[:n]):
        name, read_seq = fastq_reads[idx]
        seq1 = read_seq[:read_len]
        seq2, mutations = mutate_sequence(seq1, snp_rate, ins_rate, del_rate, error_rate)
        pairs.append({
            "index": i, "source": "real_fastq_synthetic_ref",
            "read_name": name,
            "ref_pos": -1,
            "seq1": seq1,
            "seq2": seq2,
            "mutations": mutations,
        })
    return pairs


def _build_pairs_bam(bam_path, fasta_seqs, num_pairs, read_len, min_mapq=20):
    """Build pairs from a mapped BAM: seq1=real read, seq2=reference window at true mapping site."""
    if not _PYSAM_AVAILABLE:
        raise RuntimeError(
            "pysam is required for --bam mode. Install with: pip install pysam")
    pairs = []
    bam = pysam.AlignmentFile(bam_path, "rb")
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
        seq1 = _clean_bases(raw_seq[:read_len])   # real Illumina read
        seq2 = ref[pos : pos + read_len]           # reference window at true mapping site
        pairs.append({
            "index":     len(pairs),
            "source":    "bam_mapped",
            "read_name": read.query_name,
            "ref_pos":   pos,
            "seq1":      seq1,
            "seq2":      seq2,
            "mutations": {},
        })
        if len(pairs) >= num_pairs:
            break
    bam.close()
    print(f"  BAM: {len(pairs)} pairs selected, {skipped} reads skipped "
          f"(unmapped/secondary/low-mapq/short)")
    if len(pairs) < num_pairs:
        print(f"  Warning: only {len(pairs)} usable reads found (requested {num_pairs})")
    return pairs


def _compute_scores(pairs, algorithms, bitap_max_k: int = 5):
    print(f"Computing ground truth for {len(pairs)} pairs...")
    for algo in algorithms:
        print(f"  {algo}...", end=" ", flush=True)
        for pair in pairs:
            s1, s2 = pair["seq1"], pair["seq2"]
            if   algo == "sw":        pair["score_sw"]        = smith_waterman(s1, s2)
            elif algo == "nw":        pair["score_nw"]        = needleman_wunsch(s1, s2)
            elif algo == "banded_sw": pair["score_banded_sw"] = banded_smith_waterman(s1, s2)
            elif algo == "myers":     pair["score_myers"]     = (myers_edit_distance(s1, s2)
                                                                  if len(s1) <= 32
                                                                  else edit_distance_dp(s1, s2))
            elif algo == "bp_lcs":    pair["score_bp_lcs"]    = (bitparallel_lcs(s1, s2)
                                                                  if len(s1) <= 32 else -1)
            elif algo == "wfa":       pair["score_wfa"]       = wfa_align(s1, s2)
            elif algo == "bitap":     pair["score_bitap"]     = bitap_edit_distance(s1, s2, bitap_max_k)
        print("done")


def generate_dataset(preset: TestPreset, out_dir: str,
                     algorithms: Optional[List[str]] = None, seed: int = 42,
                     fastq_path: Optional[str] = None,
                     fasta_path: Optional[str] = None,
                     bam_path: Optional[str] = None,
                     min_read_len: int = 0,
                     max_read_len: int = 0,
                     min_mapq: int = 20,
                     bitap_max_k: int = 5):
    random.seed(seed)
    os.makedirs(out_dir, exist_ok=True)

    if algorithms is None:
        algorithms = ["sw", "nw", "banded_sw", "myers", "bp_lcs", "wfa", "bitap"]

    if bam_path and fasta_path:
        print("Mode: BAM-MAPPED — real reads at true mapping positions in reference")
        fasta_seqs = read_fasta(fasta_path)
        pairs = _build_pairs_bam(bam_path, fasta_seqs,
                                 preset.num_pairs, preset.read_len,
                                 min_mapq=min_mapq)
        data_source = f"bam={bam_path}, fasta={fasta_path}"

    elif fastq_path and fasta_path:
        print("Mode: REAL — FASTQ reads + FASTA reference")
        fasta_seqs  = read_fasta(fasta_path)
        fastq_reads = read_fastq(fastq_path,
                                 max_reads = preset.num_pairs * 2,
                                 min_len   = min_read_len or 20,
                                 max_len   = max_read_len or preset.read_len)
        pairs = _build_pairs_real(fastq_reads, fasta_seqs,
                                  preset.num_pairs, preset.read_len)
        data_source = f"fastq={fastq_path}, fasta={fasta_path}"

    elif fastq_path:
        print("Mode: REAL reads + SYNTHETIC reference (mutated read as seq2)")
        fastq_reads = read_fastq(fastq_path,
                                 max_reads = preset.num_pairs * 2,
                                 min_len   = min_read_len or 20,
                                 max_len   = max_read_len or preset.read_len)
        pairs = _build_pairs_fastq_only(fastq_reads, preset.num_pairs,
                                        preset.read_len,
                                        preset.snp_rate, preset.ins_rate,
                                        preset.del_rate, preset.error_rate)
        data_source = f"fastq={fastq_path}, ref=synthetic_mutation"

    else:
        print("Mode: SYNTHETIC (random reference + simulated mutations)")
        reference = random_dna(preset.ref_len)
        pairs = _build_pairs_synthetic(preset, reference)
        with open(os.path.join(out_dir, "reference.txt"), "w") as f:
            f.write(reference + "\n")
        data_source = "synthetic"

    _compute_scores(pairs, algorithms, bitap_max_k)

    metadata = {
        "preset":      preset.name,
        "read_len":    preset.read_len,
        "num_pairs":   len(pairs),
        "algorithms":  algorithms,
        "seed":        seed,
        "data_source": data_source,
        "encoding":    "2bit_per_base_16_per_int32",
        "scoring":     "BLOSUM62_DNA",
        "gap_open":    -5,
        "gap_extend":  -1,
    }
    with open(os.path.join(out_dir, "metadata.json"), "w") as f:
        json.dump(metadata, f, indent=2)

    with open(os.path.join(out_dir, "pairs.txt"), "w") as f:
        for pair in pairs:
            info = (f"Pair {pair['index']}: source={pair['source']} "
                    f"len1={len(pair['seq1'])} len2={len(pair['seq2'])}")
            if "read_name" in pair:
                info += f" read={pair['read_name']}"
            if pair.get("ref_pos", -1) >= 0:
                info += f" ref_pos={pair['ref_pos']}"
            f.write(info + "\n")
            f.write(f"  seq1: {pair['seq1']}\n")
            f.write(f"  seq2: {pair['seq2']}\n")
            for algo in algorithms:
                key = f"score_{algo}"
                if key in pair:
                    f.write(f"  {algo}: {pair[key]}\n")
            f.write("\n")

    _write_plio_dp(pairs, algorithms, out_dir, bitap_max_k)

    print(f"\nDataset written to {out_dir}/")
    print(f"  {len(pairs)} pairs, source={data_source}")
    print(f"  Scoring: BLOSUM62_DNA (match=+2, transition=-1, transversion=-2)")
    print(f"  Encoding: 2-bit per base, 16 bases per int32")
    return pairs


BITAP_MAX_PLEN = 128  # BV_WORDS * 32; bitap_kernel.cc PM/R arrays are [NUM_BASES/MAX_K+1][4]

def _write_plio_dp(pairs, algorithms, out_dir, bitap_max_k: int = 5):
    if "bitap" in algorithms:
        for pair in pairs:
            if len(pair["seq1"]) > BITAP_MAX_PLEN:
                print(f"  WARNING: pair {pair['index']} seq1 len={len(pair['seq1'])} > "
                      f"BITAP_MAX_PLEN={BITAP_MAX_PLEN}; truncating to avoid OOB in kernel")
                pair["seq1"] = pair["seq1"][:BITAP_MAX_PLEN]
                pair["score_bitap"] = bitap_edit_distance(pair["seq1"], pair["seq2"], bitap_max_k)

    seq1_words = [len(pairs)]
    seq2_words = []
    for pair in pairs:
        s1_data, s2_data = format_plio_pair(pair["seq1"], pair["seq2"])
        seq1_words.extend(s1_data)
        seq2_words.extend(s2_data)

    for fname, words in [("seq1.txt", seq1_words), ("seq2.txt", seq2_words)]:
        with open(os.path.join(out_dir, fname), "w") as f:
            for w in words:
                if w >= (1 << 31): w -= (1 << 32)
                f.write(f"{w}\n")
            f.write("0\n")  # PLIO padding (aiesimulator TLAST workaround)

    for algo in algorithms:
        key = f"score_{algo}"
        with open(os.path.join(out_dir, f"out_score.txt"), "w") as f:
            for pair in pairs:
                f.write(f"{pair.get(key, -1)}\n")


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate test data for genomic alignment AIE kernels")
    parser.add_argument("--preset", choices=PRESETS.keys(),
                        help="Predefined configuration (synthetic)")
    parser.add_argument("--fastq", metavar="FILE",
                        help="Real Illumina reads (.fastq or .fastq.gz)")
    parser.add_argument("--ref-fasta", metavar="FILE",
                        help="Reference genome (.fa/.fasta or .gz). "
                             "PhiX174: https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1")
    parser.add_argument("--bam", metavar="FILE",
                        help="Aligned BAM (sorted+indexed). Pairs real reads with reference "
                             "windows at true mapping positions. Requires --ref-fasta. "
                             "Produces realistic ~97-99%% identity pairs. "
                             "Needs: pip install pysam")
    parser.add_argument("--min-mapq", type=int, default=20,
                        help="Minimum MAPQ to include a read (default 20; filters multi-mappers)")
    parser.add_argument("--min-read-len", type=int, default=0)
    parser.add_argument("--max-read-len", type=int, default=0)
    parser.add_argument("--ref-len",    type=int)
    parser.add_argument("--read-len",   type=int)
    parser.add_argument("--num-pairs",  type=int)
    parser.add_argument("--snp-rate",   type=float, default=0.02)
    parser.add_argument("--ins-rate",   type=float, default=0.0001)
    parser.add_argument("--del-rate",   type=float, default=0.0001)
    parser.add_argument("--error-rate", type=float, default=0.005)
    parser.add_argument("--seed",       type=int,   default=42)
    parser.add_argument("--bitap-max-k", type=int,  default=5,
                        help="Max edit distance k for bitap kernel (default 5; must match MAX_K in bitap_kernel.cc)")
    parser.add_argument("--algorithms", nargs="+",
                        default=["sw", "nw", "banded_sw", "myers", "bp_lcs", "wfa", "bitap"])
    parser.add_argument("--out-dir",    default="data")
    args = parser.parse_args()

    if args.preset:
        preset = PRESETS[args.preset]
        if args.ref_len:   preset.ref_len   = args.ref_len
        if args.read_len:  preset.read_len  = args.read_len
        if args.num_pairs: preset.num_pairs = args.num_pairs
    elif args.bam:
        read_len  = args.read_len  or 150
        num_pairs = args.num_pairs or 1000
        preset = TestPreset("bam_mapped", 0, read_len, num_pairs,
                            0.0, 0.0, 0.0, 0.0, "BAM-mapped real Illumina reads")
    elif args.fastq:
        read_len  = args.read_len  or (args.max_read_len or 150)
        num_pairs = args.num_pairs or 1000
        preset = TestPreset("real", 0, read_len, num_pairs,
                            args.snp_rate, args.ins_rate, args.del_rate,
                            args.error_rate, "Real Illumina data")
    elif args.ref_len and args.read_len and args.num_pairs:
        preset = TestPreset("custom", args.ref_len, args.read_len, args.num_pairs,
                            args.snp_rate, args.ins_rate, args.del_rate,
                            args.error_rate, "Custom synthetic")
    else:
        parser.error("Provide --preset, --fastq, or --ref-len/--read-len/--num-pairs")

    algos = list(args.algorithms)
    if preset.read_len > 32:
        algos = [a for a in algos if a != "bp_lcs"]
        if "bp_lcs" in args.algorithms:
            print(f"Note: Skipping bp_lcs (read_len {preset.read_len} > 32)")

    if args.bam and not args.ref_fasta:
        parser.error("--bam requires --ref-fasta")

    print(f"Preset: {preset.name}")
    generate_dataset(preset, args.out_dir, algos, args.seed,
                     fastq_path   = args.fastq,
                     fasta_path   = args.ref_fasta,
                     bam_path     = args.bam,
                     min_read_len = args.min_read_len,
                     max_read_len = args.max_read_len,
                     min_mapq     = args.min_mapq,
                     bitap_max_k  = args.bitap_max_k)


if __name__ == "__main__":
    main()
