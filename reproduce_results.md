# Reproduce Results

All commands below are run from `genomics-alignment/` unless stated otherwise:

```bash
cd genomics-alignment
```

For toolchain setup and general build details, see the top-level [README.md](README.md).

The paper's results were measured on an AMD VCK5000 (Versal VC1902) with an Intel Xeon E5-2640 v4 host, Vitis 2023.1, XRT 2.14.354 and Ubuntu 20.04.6 (paper §IV).

## Where each paper result comes from

| Paper result | Section below | Output file |
|---|---|---|
| Dataset creation (§IV-A) | [1. Data generation](#1-data-generation) | `test_data/illumina_real/*`, `aie/data/` |
| **Table II**: per-tile program-memory footprint | [3.2 AIEsim](#32-aiesim) (`kernel_prog_mem_bytes` column) | `results_aiesim_sweep.csv` |
| **Table III**: single-tile throughput (MCUPS) and scaling factor | [3.3 Board run: 1 tile](#33-board-run-1-tile-1x1); scaling = Table V peak ÷ Table III | `kernels/genomics_server/results_sweep.csv` |
| **Table IV**: cycles per DP cell, and §IV-B VLIW issue rate (≈0.14) and vector utilization | [3.2 AIEsim](#32-aiesim) (`cycles_per_cell`, `vliw_issue_rate`, `vec_util_pct` columns) | `results_aiesim_sweep.csv` |
| **Table V**: full-array throughput (GCUPS), 1 iteration and peak | [3.1 Board run: full system](#31-board-run-full-system-8x50) | `kernels/genomics_server/results_sweep.csv` |
| **Table VI**: performance per watt (GCUPS/W), and the ≈50 W figure | [3.5 Power](#35-power-8x50) | `kernels/genomics_server/results_power_sweep.csv` |
| **Table VII**: PL resource utilization (1x1, 4x25, 8x50) | [3.4 PL resource utilization](#34-pl-resource-utilization-one-algorithm-at-a-time-eg-wfa) | `kernels/genomics_server/results_pl_utilization.csv` |

## 0. Prerequisites

```bash
# Bioinformatics tools (SRA toolkit, samtools, bwa, wgsim)
conda install -c bioconda sra-tools samtools bwa wgsim
# or, on Debian/Ubuntu:
# sudo apt install sra-toolkit samtools bwa

# Python library for reading BAM files
pip install pysam

# Make the conda binaries visible (if not already on PATH)
export PATH=$PATH:$(conda info --base)/bin

# Verify
fastq-dump --version
samtools --version
bwa 2>&1 | head -3
wgsim 2>&1 | head -3
```

## 1. Data generation

> **Paper:** §IV-A (Dataset creation).

### a) Real dataset setup

> **Note:** `test_data/illumina_real/chr22.fa` and `test_data/illumina_real/NA12878_50k_R1.fastq.gz` are already committed to the repository, so Steps 1 and 2 are optional. Start at Step 3 if you use the committed files.

#### Step 1: Download Illumina reads (NA12878, first 50K reads)

```bash
# The NCBI FTP endpoint may be blocked, but HTTPS works.
# Stream only the first 50K reads (4 lines per read) instead of downloading the full file:
curl -s "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/\
NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz" \
  | zcat | head -200000 | gzip > test_data/illumina_real/NA12878_50k_R1.fastq.gz
```

#### Step 2: Download the reference genome (chr22 only)

Downloading all of GRCh38 is 3.2 GB. For AIE kernel testing, chromosome 22 (the smallest autosome, ~51 MB) is enough:

```bash
# Download chr22 from UCSC hg38
wget -O test_data/illumina_real/chr22.fa.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz

gunzip test_data/illumina_real/chr22.fa.gz   # produces chr22.fa (~51 MB uncompressed)

# Verify
head -2 test_data/illumina_real/chr22.fa
# >chr22
# NNNNNNNNNNNNNNNN...  (telomeric N-padding; the real sequence starts further in)
```

#### Step 3: Simulate 128 bp reads from the chr22 reference

> **Note:** The 128 bp reads used for the real-data (`bam`) experiments are created from the **chr22 reference** (`test_data/illumina_real/chr22.fa`) with `wgsim`, which gives 100,000 reads. In Step 4 they are aligned back to chr22 with `bwa mem`. The `NA12878_` prefix in the file names is historical.

```bash
wgsim -1 128 -2 128 -N 100000 -r 0.001 -R 0.10 -e 0.005 -S 42 \
    test_data/illumina_real/chr22.fa \
    test_data/illumina_real/NA12878_128bp_R1.fastq \
    test_data/illumina_real/NA12878_128bp_R2.fastq
```

Only `R1` is used in the following steps.

#### Step 4: Align the reads and create the BAM file

```bash
bwa index test_data/illumina_real/chr22.fa

bwa mem test_data/illumina_real/chr22.fa \
    test_data/illumina_real/NA12878_128bp_R1.fastq \
    | samtools sort -o test_data/illumina_real/NA12878_chr22_128bp.bam
samtools index test_data/illumina_real/NA12878_chr22_128bp.bam
```

#### Step 5: Generate AIE test data

The sweep scripts in Section 3 regenerate their own data automatically. Run the commands below only to build or simulate a single configuration by hand.

**Single-tile kernels** (`kernels/<algo>/aie/data/`). The example uses bitap with `MAX_K=5`, matching the hardware default (`#define MAX_K 5` in `kernels/bitap/aie/src/bitap_kernel.cc`), and the AIEsim sweep settings (16 pairs, 32 bp):

```bash
python3 scripts/generate_test_data.py --algorithms bitap --out-dir kernels/bitap/aie/data \
    --num-pairs 16 --ref-len 32 --read-len 32 --bitap-max-k 5
```

The graph reads these files from `kernels/<algo>/aie/data/`:

- `seq1.txt`: PLIO input stream 1 (reads)
- `seq2.txt`: PLIO input stream 2 (reference windows)
- `out_score.txt`: ground-truth scores, compared against the simulator output

**Inter-sequence SWG** (`kernels/swg_interseq/aie/data/`):

```bash
python3 scripts/gen_interseq_data.py \
    --bam test_data/illumina_real/NA12878_chr22_128bp.bam \
    --ref-fasta test_data/illumina_real/chr22.fa \
    --seq-len 32 --num-batches 2 --batch-size 16 --min-mapq 20 --seed 42 \
    --out-dir kernels/swg_interseq/aie/data
```

**genomics_server array** (`kernels/genomics_server/aie/data/`). The arguments are `<out_dir> <rows> <cols> <refs_per_row> <num_iters> <buf_size> [seed]`, and the rows/cols must match the build:

```bash
python3 kernels/genomics_server/scripts/gen_multi_data.py \
    kernels/genomics_server/aie/data 1 1 4 2 2048 42 \
    --kernel sw --read-len 128 \
    --bam test_data/illumina_real/NA12878_chr22_128bp.bam \
    --ref-fasta test_data/illumina_real/chr22.fa
# --kernel: wfa | sw | bitap   (sw is used for striped_swg, interseq_swg and scalar_swg)
```

## 2. Build kernels

The xclbin for each algorithm is built from `kernels/genomics_server/`. For the full-system results, build all five algorithms at **8x50** (and again at **1x1** for the single-tile board run).

> **Expected time:** each 8x50 hardware build (link + package) takes several hours. The 1x1 builds are much faster.

### a) Set the array size

The array size must be the same in four places. For 8x50 (the default in the repository), nothing needs to change. For 1x1:

| Location | Change for 1x1 |
|---|---|
| `kernels/genomics_server/pl/gen_servers.py` | `NUM_ROWS = 1`, `NUM_COLS = 1`, then run `cd kernels/genomics_server/pl && python3 gen_servers.py` |
| `kernels/genomics_server/hw_link/config.cfg` | replace with `config_1.cfg` (save the 8x50 version first; `config_400.cfg` is a copy of the 8x50 config) |
| `kernels/genomics_server/host/host_genomics_server.cpp` | `#define NUM_ROWS 1`, `#define NUM_COLS 1`, `#define REF_PER_ROW 4` |
| `make` command line | `ROWS=1 COLS=1` (see below) |

For 8x50, the host values are `NUM_ROWS 8`, `NUM_COLS 50` and `REF_PER_ROW 200`.

### b) Build each algorithm

```bash
cd kernels/genomics_server

# Repeat for each KERNEL in: wfa  bitap  striped_swg  interseq_swg  scalar_swg
make clean                              # required: kernel_config.h is not regenerated otherwise
make all KERNEL=wfa ROWS=8 COLS=50      # use ROWS=1 COLS=1 for the 1-tile build
```

You don't need to edit `KERNEL ?=` in `aie/Makefile`. The command-line value is passed down to the AIE sub-make.

### c) Name the xclbins

The sweep scripts look for `genomics-alignment/<algo>.xclbin`. Copy each build there right after it finishes, because the next `make clean` deletes it:

| `KERNEL=` | Copy `build.hw/genomics_server.xclbin` to |
|---|---|
| `wfa` | `genomics-alignment/genomics_server_wfa.xclbin` |
| `bitap` (`MAX_K=5`) | `genomics-alignment/genomics_server_bitap_k5.xclbin` |
| `striped_swg` | `genomics-alignment/genomics_server_striped_swg.xclbin` |
| `interseq_swg` | `genomics-alignment/genomics_server_interseq.xclbin` |
| `scalar_swg` | `genomics-alignment/genomics_server_scalar_swg.xclbin` |

For example:

```bash
cp build.hw/genomics_server.xclbin ../../genomics_server_wfa.xclbin
```

Keep the 8x50 and 1x1 xclbins in separate directories, because they share the same names. Copy the set you want to run into `genomics-alignment/` before each sweep.

## 3. Reproduce results

Both sweep scripts work out their paths from their own location, so you can launch them from any directory. Before running, check `BDF=0000:82:00.1` (in `run_power_sweep.sh`, and the default in the host) against `xbutil examine`.

### 3.1 Board run: full system (8x50)

> **Paper:** Table V (the `num_iters=1` rows give the "1 iter" columns, and the saturated rows at 100–1000 iterations give "Peak"); §IV-C.

1. Build and place the 8x50 xclbins (Section 2), and make sure the host has `NUM_ROWS 8`, `NUM_COLS 50` and `REF_PER_ROW 200`.
2. In `scripts/run_sweep.sh`, set:

   ```bash
   NUM_ROWS=8
   NUM_COLS=50
   REF_PER_ROW=200
   ```

3. Run:

   ```bash
   bash scripts/run_sweep.sh
   ```

The script patches `#define NUM_ITERS` in the host source for each value in `ITER_LIST` (1, 10, 100, 1000), rebuilds the host, generates synthetic and BAM data for every algorithm, and runs on the board. When it finishes, it resets `NUM_ITERS` to 1.

Results are written to `kernels/genomics_server/results_sweep.csv`, and per-run logs go to `kernels/genomics_server/logs_sweep/`.

### 3.2 AIEsim

> **Paper:** Table II (`kernel_prog_mem_bytes`), Table IV (`cycles_per_cell`) and the §IV-B profile figures (`vliw_issue_rate`, `vec_util_pct`), using 32 bp sequences and 16 pairs on both datasets.

```bash
python3 scripts/aiesim_sweep.py
```

This builds and simulates each single-tile kernel in `kernels/<algo>/aie/` with `aiesimulator`. Options: `--dry-run`, `--algo bitap wfa`, `--output <file>.csv`.

The sweep values are set under "Sweep configuration" at the top of `scripts/aiesim_sweep.py`. The trailing comment on each line lists the other values you can use:

```python
DATA_TYPES  = ["synthetic"]   # options: ["synthetic", "bam"]
NUM_PAIRS   = [16]            # options: [4, 8, 16]
READ_LENS   = [32]            # options: [64, 128]
```

For the paper's numbers, use `DATA_TYPES = ["synthetic", "bam"]`, `NUM_PAIRS = [16]` and `READ_LENS = [32]`. Larger `NUM_PAIRS` and `READ_LENS` values make the simulation much slower.

Results are written to `results_aiesim_sweep.csv`.

### 3.3 Board run: 1 tile (1x1)

> **Paper:** Table III (single-tile MCUPS). The scaling factor is the Table V peak divided by this rate; §IV-B.

1. Build and place the 1x1 xclbins (Section 2), and make sure the host has `NUM_ROWS 1`, `NUM_COLS 1` and `REF_PER_ROW 4`.
2. In `scripts/run_sweep.sh`, set:

   ```bash
   NUM_ROWS=1
   NUM_COLS=1
   REF_PER_ROW=4
   ```

   `REF_PER_ROW=4` puts 4 pairs on the single tile. This matches the full array's load (200 refs per row ÷ 50 columns = 4 pairs per tile), so the Table III rates can be compared directly with Table V.

3. Run:

   ```bash
   bash scripts/run_sweep.sh
   ```

Results are written to `kernels/genomics_server/results_sweep.csv`.

### 3.4 PL resource utilization (one algorithm at a time, e.g. WFA)

> **Paper:** Table VII (WFA at 1x1, 4x25 and 8x50); §IV-D.

The script reads `kernels/genomics_server/aie/build.hw/kernel_config.h` and the reports of the **most recent** hardware build. It checks that the build matches `EXPECTED_ROWS` / `EXPECTED_COLS`, set near the top of `scripts/gen_pl_utilization_report.py`. Each run adds or updates one row in the CSV, so you can build up the table one build at a time.

#### a) 1x1

1. Build with `make all KERNEL=wfa ROWS=1 COLS=1` (Section 2).
2. In `scripts/gen_pl_utilization_report.py`, set `EXPECTED_ROWS = 1` and `EXPECTED_COLS = 1`.
3. Run:

   ```bash
   python3 scripts/gen_pl_utilization_report.py --algo genomics_server_wfa
   ```

#### b) 4x25

Repeat with `ROWS=4 COLS=25` and `EXPECTED_ROWS = 4`, `EXPECTED_COLS = 25`. There is no 4x25 link config in the repository, so first create `hw_link/config.cfg` with `in_0`..`in_3` and `out_0`..`out_24`, and set `NUM_ROWS = 4`, `NUM_COLS = 25` in `gen_servers.py`.

#### c) 8x50

Repeat with `ROWS=8 COLS=50` and `EXPECTED_ROWS = 8`, `EXPECTED_COLS = 50`.

For bitap, always pass `--algo genomics_server_bitap_k5` explicitly. The script can't tell the k5 and k100 builds apart from `kernel_config.h`.

Results are written to `kernels/genomics_server/results_pl_utilization.csv` (use `--list` to print it).

### 3.5 Power (8x50)

> **Paper:** Table VI (GCUPS/W over 100 batches) and the ≈50 W board power; §IV-C. Power is sampled with `xbutil` over a 120 s window.

1. Use the 8x50 xclbins and host settings (`NUM_ROWS 8`, `NUM_COLS 50`, `REF_PER_ROW 200`).
2. In `kernels/genomics_server/host/host_genomics_server.cpp`, set `#define TOTAL_ITER 100`. The script sets `NUM_ITERS 100` itself, but it does **not** change `TOTAL_ITER`.
3. Run:

   ```bash
   bash scripts/run_power_sweep.sh
   ```

4. When you're done, set `TOTAL_ITER` back to `1` for the throughput runs.

Results are written to `kernels/genomics_server/results_power_sweep.csv`, and logs go to `kernels/genomics_server/logs_power_sweep/`.
