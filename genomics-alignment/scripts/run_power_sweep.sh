#!/bin/bash
# Power sweep: data_type × algo/kernel
# Rebuilds the host with NUM_ITERS=100, TOTAL_ITER=100 before sweeping, so it
# never runs against a stale binary left over from run_sweep.sh (which sweeps
# NUM_ITERS up to 1000) or a manual edit.
#
# Usage: bash scripts/run_power_sweep.sh
# Run from: project root (genomics-alignment/)

set -uo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
HOST_DIR=$ROOT/kernels/genomics_server/host
HOST_CPP=$HOST_DIR/host_genomics_server.cpp
DATA_DIR=$ROOT/kernels/genomics_server/aie/data
SCRIPT=$ROOT/kernels/genomics_server/scripts/gen_multi_data.py
POWER_SCRIPT=$ROOT/scripts/monitor_power.py
CSV=$ROOT/kernels/genomics_server/results_power_sweep.csv
LOG_DIR=$ROOT/kernels/genomics_server/logs_power_sweep
BDF=0000:82:00.1
FASTQ=$ROOT/test_data/illumina_real/NA12878_50k_R1.fastq.gz
REF_FASTA=$ROOT/test_data/illumina_real/chr22.fa
BAM=$ROOT/test_data/illumina_real/NA12878_chr22_128bp.bam
MIN_MAPQ=20

mkdir -p "$LOG_DIR"

# ── rebuild host with the NUM_ITERS this script actually generates data for ──
echo "Rebuilding host with NUM_ITERS=100 ..."
sed -i "s/#define NUM_ITERS[[:space:]].*/#define NUM_ITERS   100/" "$HOST_CPP"
cd "$HOST_DIR"
make clean > /dev/null 2>&1
if ! make all > "$LOG_DIR/build_host.log" 2>&1; then
    echo "BUILD FAILED — see $LOG_DIR/build_host.log"
    exit 1
fi
echo "  Build OK (log: $LOG_DIR/build_host.log)"

# xclbin name → --kernel flag for gen_multi_data.py
declare -A ALGO_KERNEL=(
    [genomics_server_bitap_k5]=bitap
    [genomics_server_interseq]=sw
    [genomics_server_scalar_swg]=sw
    [genomics_server_striped_swg]=sw
    [genomics_server_wfa]=wfa
)
ALGO_ORDER=(
    genomics_server_bitap_k5
    genomics_server_interseq
    genomics_server_scalar_swg
    genomics_server_striped_swg
    genomics_server_wfa
)
DATA_TYPES=(synthetic bam)

echo "data_type,algo,kernel,pass_fail,avg_w,min_w,median_w,max_w,gcups,gcups_per_watt" > "$CSV"

# Track which (dt, kernel) combinations have already had data generated.
declare -A datagen_done

for dt in "${DATA_TYPES[@]}"; do
    echo "================================================================"
    echo "data_type: $dt"

    for algo in "${ALGO_ORDER[@]}"; do
        kernel="${ALGO_KERNEL[$algo]}"
        printf "\n  algo=%-38s  kernel=%s\n" "$algo" "$kernel"

        XCLBIN=$ROOT/${algo}.xclbin
        if [ ! -f "$XCLBIN" ]; then
            echo "    NO_XCLBIN — skipping"
            echo "$dt,$algo,$kernel,NO_XCLBIN,N/A,N/A,N/A,N/A,N/A,N/A" >> "$CSV"
            continue
        fi

        # ── data generation (once per dt+kernel combination) ─────────────────
        datagen_key="${dt}_${kernel}"
        if [ -z "${datagen_done[$datagen_key]+x}" ]; then
            cd "$ROOT"
            datagen_log="$LOG_DIR/datagen_${dt}_${kernel}.log"

            MAX_K_ARG=""
            [ "$kernel" = "bitap" ] && MAX_K_ARG="--bitap-max-k 5"

            if [ "$dt" = "synthetic" ]; then
                python3 "$SCRIPT" "$DATA_DIR" 8 50 200 100 2048 42 \
                    --kernel "$kernel" --read-len 128 $MAX_K_ARG \
                    > "$datagen_log" 2>&1
            elif [ "$dt" = "bam" ]; then
                python3 "$SCRIPT" "$DATA_DIR" 8 50 200 100 2048 42 \
                    --kernel "$kernel" --read-len 128 $MAX_K_ARG \
                    --bam "$BAM" --ref-fasta "$REF_FASTA" --min-mapq "$MIN_MAPQ" \
                    > "$datagen_log" 2>&1
            else
                python3 "$SCRIPT" "$DATA_DIR" 8 50 200 100 2048 42 \
                    --kernel "$kernel" --read-len 128 $MAX_K_ARG \
                    --fastq "$FASTQ" --ref-fasta "$REF_FASTA" \
                    > "$datagen_log" 2>&1
            fi

            if [ $? -ne 0 ]; then
                echo "    DATAGEN_FAIL — see $datagen_log"
                echo "$dt,$algo,$kernel,DATAGEN_FAIL,N/A,N/A,N/A,N/A,N/A,N/A" >> "$CSV"
                continue
            fi
            datagen_done[$datagen_key]=1
        else
            echo "    Data reused (kernel=$kernel already generated for $dt)"
        fi

        # ── start power monitor in background ────────────────────────────────
        # Large -t so it outlives the hw run; we kill it when hw finishes.
        power_log="$LOG_DIR/power_${dt}_${algo}.log"
        python3 "$POWER_SCRIPT" -d "$BDF" -t 120 > "$power_log" 2>&1 &
        POWER_PID=$!
        echo "    Power monitor started (pid=$POWER_PID)"

        # ── hw run ───────────────────────────────────────────────────────────
        cd "$HOST_DIR"
        run_log="$LOG_DIR/run_${dt}_${algo}.log"
        ./genomics_server "$XCLBIN" "$DATA_DIR" > "$run_log" 2>&1
        RUN_RC=$?

        # ── stop power monitor ────────────────────────────────────────────────
        kill "$POWER_PID" 2>/dev/null
        wait "$POWER_PID" 2>/dev/null || true

        # ── PASS/FAIL ─────────────────────────────────────────────────────────
        if [ $RUN_RC -ne 0 ]; then
            pf=RUN_ERROR
        elif grep -q "^PASS$" "$run_log"; then
            pf=PASS
        else
            pf=FAIL
        fi

        # ── parse power stats from individual readings ────────────────────────
        # monitor_power.py prints "  Power : 45 Watts" for each sample.
        # We killed it before the summary, so compute stats from raw lines.
        read -r avg_w min_w median_w max_w < <(python3 - <<PYEOF
import re, statistics, sys
vals = []
with open("$power_log") as f:
    for line in f:
        m = re.search(r'Power\s*:\s*(\d+)\s*Watts', line)
        if m:
            vals.append(int(m.group(1)))
if vals:
    print(f"{statistics.mean(vals):.2f}", min(vals),
          statistics.median(vals), max(vals))
else:
    print("N/A N/A N/A N/A")
PYEOF
)

        # ── extract throughput ─────────────────────────────────────────────────
        # Matches: "Throughput:       8.42 GCUPS  [wall clock]"
        gcups=$(grep -oP 'Throughput:\s+\K[\d.]+(?=\s+GCUPS)' "$run_log" | head -1)
        [ -z "$gcups" ] && gcups=N/A

        # ── performance per watt (GCUPS / avg_w) ────────────────────────────────
        if [ "$gcups" != "N/A" ] && [ "$avg_w" != "N/A" ]; then
            perf_per_watt=$(python3 -c "print(f'{$gcups / $avg_w:.4f}')")
        else
            perf_per_watt=N/A
        fi

        printf "    %-10s  avg=%sW  min=%sW  median=%sW  max=%sW  %sGCUPS  %sGCUPS/W\n" \
            "$pf" "$avg_w" "$min_w" "$median_w" "$max_w" "$gcups" "$perf_per_watt"
        echo "$dt,$algo,$kernel,$pf,$avg_w,$min_w,$median_w,$max_w,$gcups,$perf_per_watt" >> "$CSV"
    done
done

echo ""
echo "Power sweep complete."
echo "Results : $CSV"
echo "Logs    : $LOG_DIR/"
