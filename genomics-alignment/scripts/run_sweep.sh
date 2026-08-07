#!/bin/bash
# Sweep: data_type × num_iters × algo/kernel
# Outputs results_sweep.csv and per-run logs under logs_sweep/
#
# Usage: bash scripts/run_sweep.sh
# Run from: project root (genomics-alignment/)

set -uo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
HOST_DIR=$ROOT/kernels/genomics_server/host
HOST_CPP=$HOST_DIR/host_genomics_server.cpp
DATA_DIR=$ROOT/kernels/genomics_server/aie/data
SCRIPT=$ROOT/kernels/genomics_server/scripts/gen_multi_data.py
CSV=$ROOT/kernels/genomics_server/results_sweep.csv
LOG_DIR=$ROOT/kernels/genomics_server/logs_sweep
FASTQ=$ROOT/test_data/illumina_real/NA12878_50k_R1.fastq.gz
REF_FASTA=$ROOT/test_data/illumina_real/chr22.fa
BAM=$ROOT/test_data/illumina_real/NA12878_chr22_128bp.bam
MIN_MAPQ=20
NUM_ROWS=1
NUM_COLS=1
REF_PER_ROW=32

mkdir -p "$LOG_DIR"

# xclbin name → --kernel flag for gen_multi_data.py
declare -A ALGO_KERNEL=(
    [genomics_server_bitap_k5]=bitap
    [genomics_server_interseq]=sw
    [genomics_server_scalar_swg]=sw
    [genomics_server_striped_swg]=sw
    [genomics_server_wfa]=wfa
)
# Fixed order for CSV readability
ALGO_ORDER=(
    genomics_server_bitap_k5
    genomics_server_interseq
    genomics_server_scalar_swg
    genomics_server_striped_swg
    genomics_server_wfa
)

ITER_LIST=(1 10 100 1000) #(1 10 100 1000)
DATA_TYPES=(synthetic bam) #(synthetic bam)

echo "data_type,num_iters,algo,kernel,pass_fail,throughput_gcups" > "$CSV"

for iters in "${ITER_LIST[@]}"; do
    echo "================================================================"
    echo "Rebuilding host with NUM_ITERS=$iters ..."

    # Patch #define NUM_ITERS in-place (matches the line regardless of whitespace)
    sed -i "s/#define NUM_ITERS[[:space:]].*/#define NUM_ITERS   $iters/" "$HOST_CPP"

    cd "$HOST_DIR"
    make clean > /dev/null 2>&1
    if ! make all > "$LOG_DIR/build_iters${iters}.log" 2>&1; then
        echo "  BUILD FAILED — see $LOG_DIR/build_iters${iters}.log"
        for dt in "${DATA_TYPES[@]}"; do
            for algo in "${ALGO_ORDER[@]}"; do
                kernel="${ALGO_KERNEL[$algo]}"
                echo "$dt,$iters,$algo,$kernel,BUILD_FAIL,N/A" >> "$CSV"
            done
        done
        continue
    fi
    echo "  Build OK (log: $LOG_DIR/build_iters${iters}.log)"

    for dt in "${DATA_TYPES[@]}"; do
        echo "  --- data_type=$dt  num_iters=$iters ---"

        for algo in "${ALGO_ORDER[@]}"; do
            kernel="${ALGO_KERNEL[$algo]}"
            printf "    %-38s kernel=%-12s " "$algo" "$kernel"

            # ── data generation ──────────────────────────────────────────────
            cd "$ROOT"
            datagen_log="$LOG_DIR/datagen_${dt}_i${iters}_${algo}.log"

            # --bitap-max-k only matters when kernel=bitap; harmless for others
            MAX_K_ARG=""
            if [ "$kernel" = "bitap" ]; then
                MAX_K_ARG="--bitap-max-k 5"
            fi

            if [ "$dt" = "synthetic" ]; then
                python3 "$SCRIPT" "$DATA_DIR" "$NUM_ROWS" "$NUM_COLS" "$REF_PER_ROW" "$iters" 2048 42 \
                    --kernel "$kernel" --read-len 128 $MAX_K_ARG \
                    > "$datagen_log" 2>&1
            elif [ "$dt" = "bam" ]; then
                python3 "$SCRIPT" "$DATA_DIR" "$NUM_ROWS" "$NUM_COLS" "$REF_PER_ROW" "$iters" 2048 42 \
                    --kernel "$kernel" --read-len 128 $MAX_K_ARG \
                    --bam "$BAM" --ref-fasta "$REF_FASTA" --min-mapq "$MIN_MAPQ" \
                    > "$datagen_log" 2>&1
            else
                python3 "$SCRIPT" "$DATA_DIR" "$NUM_ROWS" "$NUM_COLS" "$REF_PER_ROW" "$iters" 2048 42 \
                    --kernel "$kernel" --read-len 128 $MAX_K_ARG \
                    --fastq "$FASTQ" --ref-fasta "$REF_FASTA" \
                    > "$datagen_log" 2>&1
            fi

            if [ $? -ne 0 ]; then
                echo "DATAGEN_FAIL"
                echo "$dt,$iters,$algo,$kernel,DATAGEN_FAIL,N/A" >> "$CSV"
                continue
            fi

            # ── hw run ───────────────────────────────────────────────────────
            XCLBIN=$ROOT/${algo}.xclbin
            if [ ! -f "$XCLBIN" ]; then
                echo "NO_XCLBIN"
                echo "$dt,$iters,$algo,$kernel,NO_XCLBIN,N/A" >> "$CSV"
                continue
            fi

            cd "$HOST_DIR"
            run_log="$LOG_DIR/run_${dt}_i${iters}_${algo}.log"
            ./genomics_server "$XCLBIN" "$DATA_DIR" > "$run_log" 2>&1

            # ── extract PASS/FAIL ────────────────────────────────────────────
            if grep -q "^PASS$" "$run_log"; then
                pf=PASS
            else
                pf=FAIL
            fi

            # ── extract throughput ───────────────────────────────────────────
            # Matches: "Throughput:       8.42 GCUPS  [wall clock]"
            gcups=$(grep -oP 'Throughput:\s+\K[\d.]+(?=\s+GCUPS)' "$run_log" | head -1)
            [ -z "$gcups" ] && gcups=N/A

            printf "%-6s  %s GCUPS\n" "$pf" "$gcups"
            echo "$dt,$iters,$algo,$kernel,$pf,$gcups" >> "$CSV"
        done
    done
done

# Restore NUM_ITERS to 1 (default)
sed -i "s/#define NUM_ITERS[[:space:]].*/#define NUM_ITERS   1/" "$HOST_CPP"

echo ""
echo "Sweep complete."
echo "Results : $CSV"
echo "Logs    : $LOG_DIR/"
