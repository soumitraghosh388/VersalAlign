#!/usr/bin/env python3
"""Generate test data for inter-sequence SIMD SWG kernel.

Batches BATCH_SIZE pairs together (all same len1 x len2).
Output format matches swg_interseq_kernel.cc protocol (2-bit encoding).

Synthetic mode (default):
    python gen_interseq_data.py --seq-len 32 --num-batches 2

Real Illumina mode (FASTQ only — seq2 is a mutated version of the read):
    python gen_interseq_data.py --fastq reads.fastq.gz --seq-len 150

Real Illumina mode (FASTQ + reference FASTA):
    python gen_interseq_data.py --fastq reads.fastq.gz --ref-fasta chr22.fa --seq-len 150

Real BAM mode (mapped reads + reference at true mapping positions):
    python gen_interseq_data.py --bam reads.bam --ref-fasta chr22.fa --seq-len 150
"""

import argparse
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', 'scripts'))
from generate_test_data import (
    smith_waterman, BASE_TO_INT, BASES,
    read_fastq, read_fasta, fasta_window,
    _clean_bases, mutate_sequence, pack_bases_to_int32,
    _build_pairs_bam,
)


def signed32(v):
    if v >= (1 << 31):
        v -= (1 << 32)
    return v


def build_pairs_synthetic(seq_len, batch_size, num_batches, snp_rate, ins_rate, del_rate, error_rate):
    all_pairs = []
    for _ in range(num_batches * batch_size):
        ref = ''.join(random.choice(BASES) for _ in range(seq_len))
        read, _ = mutate_sequence(ref, snp_rate, ins_rate, del_rate, error_rate)
        read = read[:seq_len].ljust(seq_len, random.choice(BASES))
        all_pairs.append((read, ref))
    return all_pairs


def build_pairs_fastq_only(fastq_path, seq_len, batch_size, num_batches,
                            snp_rate, ins_rate, del_rate, error_rate):
    needed = batch_size * num_batches
    reads = read_fastq(fastq_path, max_reads=needed * 4, min_len=1)
    if len(reads) < needed:
        print(f"  Warning: only {len(reads)} reads available, needed {needed}")
    random.shuffle(reads)
    all_pairs = []
    for _, seq in reads[:needed]:
        seq1 = _clean_bases(seq)[:seq_len].ljust(seq_len, random.choice(BASES))
        seq2, _ = mutate_sequence(seq1, snp_rate, ins_rate, del_rate, error_rate)
        seq2 = seq2[:seq_len].ljust(seq_len, random.choice(BASES))
        all_pairs.append((seq1, seq2))
    while len(all_pairs) < needed:
        ref = ''.join(random.choice(BASES) for _ in range(seq_len))
        read, _ = mutate_sequence(ref, snp_rate, ins_rate, del_rate, error_rate)
        all_pairs.append((read[:seq_len], ref))
    return all_pairs


def build_pairs_fastq_fasta(fastq_path, fasta_path, seq_len, batch_size, num_batches,
                             snp_rate, ins_rate, del_rate, error_rate):
    needed = batch_size * num_batches
    fasta_seqs = read_fasta(fasta_path)
    reads = read_fastq(fastq_path, max_reads=needed * 4, min_len=1)
    if len(reads) < needed:
        print(f"  Warning: only {len(reads)} reads available, needed {needed}")
    random.shuffle(reads)
    all_pairs = []
    for _, seq in reads[:needed]:
        seq1 = _clean_bases(seq)[:seq_len].ljust(seq_len, random.choice(BASES))
        _, ref_seq, _ = fasta_window(fasta_seqs, seq_len)
        all_pairs.append((seq1, ref_seq))
    while len(all_pairs) < needed:
        _, ref_seq, _ = fasta_window(fasta_seqs, seq_len)
        read, _ = mutate_sequence(ref_seq, snp_rate, ins_rate, del_rate, error_rate)
        all_pairs.append((read[:seq_len], ref_seq))
    return all_pairs


def build_pairs_bam(bam_path, fasta_path, seq_len, batch_size, num_batches, min_mapq=20):
    needed = batch_size * num_batches
    fasta_seqs = read_fasta(fasta_path)
    bam_pairs = _build_pairs_bam(bam_path, fasta_seqs, needed, seq_len, min_mapq=min_mapq)
    return [(p["seq1"], p["seq2"]) for p in bam_pairs]


def write_outputs(all_pairs, batch_size, num_batches, seq_len, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    all_scores = [smith_waterman(s1, s2) for s1, s2 in all_pairs]

    seq1_words = [num_batches]
    seq2_words = []
    lens_word = (seq_len << 16) | seq_len

    for b in range(num_batches):
        batch = all_pairs[b * batch_size:(b + 1) * batch_size]

        seq1_words.append(lens_word)
        for s1, _ in batch:
            seq1_words.extend(pack_bases_to_int32(s1))

        seq2_words.append(lens_word)
        for _, s2 in batch:
            seq2_words.extend(pack_bases_to_int32(s2))

    with open(os.path.join(out_dir, 'seq1.txt'), 'w') as f:
        for w in seq1_words:
            f.write(f'{signed32(w)}\n')
        f.write('0\n')

    with open(os.path.join(out_dir, 'seq2.txt'), 'w') as f:
        for w in seq2_words:
            f.write(f'{signed32(w)}\n')
        f.write('0\n')

    # Golden scores: 8 int32 per batch (2 int16 packed per int32, matches kernel output)
    with open(os.path.join(out_dir, 'out_score.txt'), 'w') as f:
        for b in range(num_batches):
            batch_scores = all_scores[b * batch_size:(b + 1) * batch_size]
            for p in range(0, batch_size, 2):
                s0 = batch_scores[p] & 0xFFFF
                s1 = batch_scores[p + 1] & 0xFFFF
                #f.write(f'{signed32(s0 | (s1 << 16))}\n')
                f.write(f'{signed32(s0)}\n')
                f.write(f'{signed32(s1)}\n')

    with open(os.path.join(out_dir, 'scores_detail.txt'), 'w') as f:
        for b in range(num_batches):
            f.write(f'Batch {b}:\n')
            for p in range(batch_size):
                idx = b * batch_size + p
                s1, s2 = all_pairs[idx]
                f.write(f'  Pair {p}: score={all_scores[idx]}  seq1={s1}  seq2={s2}\n')

    total = num_batches * batch_size
    words_per_seq = (seq_len + 15) // 16
    print(f'Generated {num_batches} batches × {batch_size} pairs = {total} total')
    print(f'Sequence length: {seq_len} bp  ({words_per_seq} int32/seq, 2-bit encoded)')
    print(f'seq1.txt: {len(seq1_words)+1} lines, seq2.txt: {len(seq2_words)+1} lines')
    print(f'out_score.txt: {num_batches * batch_size // 2} lines (packed int32)')
    print(f'Scores: {all_scores}')


def main():
    parser = argparse.ArgumentParser(
        description='Generate inter-sequence SWG test data (2-bit encoded, BLOSUM62)')
    parser.add_argument('--fastq',       metavar='FILE', help='Illumina reads (.fastq or .fastq.gz)')
    parser.add_argument('--ref-fasta',   metavar='FILE', help='Reference genome (.fa/.fasta or .gz)')
    parser.add_argument('--bam',         metavar='FILE', help='Mapped BAM — seq1=read, seq2=reference '
                                                               'window at true mapping position. Requires --ref-fasta.')
    parser.add_argument('--min-mapq',    type=int,   default=20,   help='Minimum MAPQ to include a read (default 20)')
    parser.add_argument('--seq-len',     type=int,   default=32,   help='Sequence length (default 32)')
    parser.add_argument('--batch-size',  type=int,   default=16,   help='Pairs per batch (default 16)')
    parser.add_argument('--num-batches', type=int,   default=2,    help='Number of batches (default 2)')
    parser.add_argument('--seed',        type=int,   default=42,   help='Random seed (default 42)')
    parser.add_argument('--snp-rate',    type=float, default=0.02, help='SNP rate for synthetic/mutated mode')
    parser.add_argument('--ins-rate',    type=float, default=0.0001,help='Insertion rate')
    parser.add_argument('--del-rate',    type=float, default=0.0001,help='Deletion rate')
    parser.add_argument('--error-rate',  type=float, default=0.005,help='Sequencing error rate')
    parser.add_argument('--out-dir',     default='.',              help='Output directory (default .)')
    args = parser.parse_args()

    random.seed(args.seed)

    if args.bam and not args.ref_fasta:
        parser.error('--bam requires --ref-fasta')

    if args.bam:
        print('Mode: BAM-MAPPED — real reads at true mapping positions in reference')
        pairs = build_pairs_bam(
            args.bam, args.ref_fasta,
            args.seq_len, args.batch_size, args.num_batches,
            min_mapq=args.min_mapq)
    elif args.fastq and args.ref_fasta:
        print('Mode: REAL — FASTQ reads paired with FASTA reference windows')
        pairs = build_pairs_fastq_fasta(
            args.fastq, args.ref_fasta,
            args.seq_len, args.batch_size, args.num_batches,
            args.snp_rate, args.ins_rate, args.del_rate, args.error_rate)
    elif args.fastq:
        print('Mode: REAL reads — seq2 is a mutated version of each read')
        pairs = build_pairs_fastq_only(
            args.fastq, args.seq_len, args.batch_size, args.num_batches,
            args.snp_rate, args.ins_rate, args.del_rate, args.error_rate)
    else:
        print('Mode: SYNTHETIC — random reference + mutations')
        pairs = build_pairs_synthetic(
            args.seq_len, args.batch_size, args.num_batches,
            args.snp_rate, args.ins_rate, args.del_rate, args.error_rate)

    write_outputs(pairs, args.batch_size, args.num_batches, args.seq_len, args.out_dir)


if __name__ == '__main__':
    main()
