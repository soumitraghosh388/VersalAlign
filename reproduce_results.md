# Reproduce Results

## 1. Data generation

### a) Real dataset setup

#### Genome ILLUMINA genome dataset download steps

```bash
# SRA Toolkit (for downloading reads from NCBI)
conda install -c bioconda sra-tools    # or
apt install sra-toolkit                 # Debian/Ubuntu

# samtools (for working with reference/BAM)
conda install -c bioconda samtools

# get conda base path
conda info --base

export PATH=$PATH:/home/soumitra/official/phd/software/anaconda3/bin

# Verify
fastq-dump --version
samtools --version
```

#### Step 2: download read

```bash
# The NCBI FTP is blocked but HTTPS works
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/\
Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz

# If still slow, stream only first 50K reads without downloading full file:
curl -s "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/\
NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz" \
  | zcat | head -200000 | gzip > NA12878_50k_R1.fastq.gz
```

#### Step 3: download reference genome (chr22 only)

Downloading all of GRCh38 is 3.2 GB. For AIE kernel testing, chromosome 22 (smallest autosome, ~51 MB) is sufficient:

```bash
# Download chr22 from UCSC hg38
wget -O chr22.fa.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz

gunzip chr22.fa.gz   # produces chr22.fa (~51MB uncompressed)

# Verify
head -2 chr22.fa
# >chr22
# NNNNNNNNNNNNNNNN...  (telomeric N-padding, real sequence starts further in)
```

#### Step 4: create reads of length 128 bp from the real reference dataset

```bash
wgsim -1 128 -2 128 -N 100000 -r 0.001 -R 0.10 -e 0.005 -S 42 \
    test_data/illumina_real/chr22.fa \
    test_data/illumina_real/NA12878_128bp_R1.fastq test_data/illumina_real/NA12878_128bp_R2.fastq
```

#### Step 5: create BAM file for read sequences matching the reference sequence

```bash
# Tools
sudo apt install bwa samtools       # or via conda

# Python library
pip install pysam
bwa index test_data/illumina_real/chr22.fa

bwa mem test_data/illumina_real/chr22.fa \
    test_data/illumina_real/NA12878_128bp_R1.fastq \
    | samtools sort -o test_data/illumina_real/NA12878_chr22_128bp.bam
samtools index test_data/illumina_real/NA12878_chr22_128bp.bam
```

#### Step 6: generate AIE test data (with BAM file)

```bash
python3 scripts/generate_test_data.py --algorithms bitap --out-dir kernels/bitap/aie/data --num-pairs 5 --ref-len 15 --read-len 15 --bitap-max-k 100
```

Output files:

- `data/na12878_chr22/seq1.txt` — PLIO input stream 1 (reads)
- `data/na12878_chr22/seq2.txt` — PLIO input stream 2 (reference windows)
- `data/na12878_chr22/out_score_sw.txt` — ground-truth SWG scores
- `data/na12878_chr22/pairs.txt` — human-readable details
- `data/na12878_chr22/metadata.json` — dataset description

```bash
python3 kernels/genomics_server/scripts/gen_multi_data.py kernels/genomics_server/aie/data 1 1 32 2 2048 42 --kernel sw --read-len 128 --bam test_data/illumina_real/NA12878_chr22_128bp.bam --ref-fasta test_data/illumina_real/chr22.fa

python3 scripts/gen_interseq_data.py --bam test_data/illumina_real/NA12878_chr22.bam --ref-fasta test_data/illumina_real/chr22.fa --seq-len 32 --num-batches 2 --out-dir kernels/swg_interseq/aie/data --min-mapq 20 --seed 42 --batch-size 16
```

## 2. Build kernels

### a) Go to the genomics server AIE makefile

```bash
cd genomics-alignment/kernels/genomics_server/aie
```

Change the kernel selection in the makefile:

```makefile
KERNEL ?= wfa
```

Then build:

```bash
cd ../..
make clean
make all
```

## 3. Results reproduce

### 1. Board run full system

- Check the host code for `num_rows=8`, `cols=50`, and `ref_per_row=200`.
- Go to the scripts directory:

```bash
cd genomics-alignment/scripts
```

- Check `run_sweep.sh` for `num_rows=8`, `cols=50`, and `ref_per_row=200`.
- Run:

```bash
./run_sweep.sh
```

- Results are written to:

```text
genomics-alignment/kernels/genomics_server/results_sweep.csv
```

### 2. AIEsim

```bash
cd genomics-alignment
python3 scripts/aiesim_sweep.py
```

- Results are written to:

```text
genomics-alignment/results_aiesim_sweep.csv
```

### 3. Board run 1 tile

- Check the host code for `num_rows=1`, `cols=1`, and `ref_per_row=32`.
- Go to the scripts directory:

```bash
cd genomics-alignment/scripts
```

- Check `run_sweep.sh` for `num_rows=1`, `cols=1`, and `ref_per_row=32`.
- Ensure that the 1-tile xclbins are present; otherwise update the path in `run_sweep.sh` near line 110.
- Run:

```bash
./run_sweep.sh
```

- Results are written to:

```text
genomics-alignment/kernels/genomics_server/results_sweep.csv
```

### 4. PL resource utilization (for one algorithm, e.g. WFA)

#### a) 1x1

- Build the design.
- In `scripts/gen_pl_utilization_report.py`, change `EXPECTED_ROWS` and `EXPECTED_COLS` to `1` and `1`.
- Run:

```bash
python3 scripts/gen_pl_utilization_report.py --algo genomics_server_wfa
```

- Results are written to:

```text
genomics-alignment/kernels/genomics_server/results_pl_utilization.csv
```

#### b) 4x25

Repeat the step above after changing `EXPECTED_ROWS` and `EXPECTED_COLS` to `4` and `25`.

#### c) 8x50

Repeat the step above after changing `EXPECTED_ROWS` and `EXPECTED_COLS` to `8` and `50`.

### 5. Power

- Check the host code for `num_rows=8`, `cols=50`, `ref_per_row=200`, `num_iter=100`, and `total_iters=100`.
- Go to the scripts directory:

```bash
cd genomics-alignment/scripts
```

- Run:

```bash
./run_power_sweep.sh
```
