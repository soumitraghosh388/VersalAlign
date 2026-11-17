#!/usr/bin/env python3
"""
aiesim_sweep.py -- AIESIM performance sweep for genomics alignment kernels.

Sweeps:  data_type × algo × num_pairs × read_len
For each combination:
  1. Generates test data via generate_test_data.py
  2. Builds AIE kernel (clean build dir + make all)
  3. Runs AIESIM  (make sim)
  4. Checks PASS / FAIL
  5. Parses profile_funct / profile_instr / AIESimulator.log / Map_Report.csv
  6. Writes all metrics to a CSV file (incremental — safe to interrupt)

Usage:
  cd /home/soumitra.ghosh/versal/vck5000/fpga_genomics_exp_git/genomics-alignment
  python3 scripts/aiesim_sweep.py                        # full sweep
  python3 scripts/aiesim_sweep.py --dry-run              # print combos, no build
  python3 scripts/aiesim_sweep.py --algo bitap wfa       # subset of algos
  python3 scripts/aiesim_sweep.py --output my_results.csv
"""

import argparse
import csv
import re
import subprocess
import time
from datetime import datetime
from pathlib import Path
from typing import Optional

# ─────────────────────────────────────────────────────────────────────────────
# Sweep configuration
# ─────────────────────────────────────────────────────────────────────────────

ROOT = Path(__file__).resolve().parent.parent   # genomics-alignment/

# algo_name -> (data_gen_algorithm, kernel_dir)
ALGO_KERNEL_MAP = {
    "genomics_server_bitap_k5":   ("bitap",       "bitap"),
    "genomics_server_interseq":     ("swg_interseq", "swg_interseq"),
    "genomics_server_scalar_swg":   ("sw", "swg_scalar"),
    "genomics_server_striped_swg":  ("sw", "sw_striped"),
    "genomics_server_wfa":          ("wfa",         "wfa"),
}

DATA_TYPES  = ["synthetic", "bam"] #["synthetic", "illumina", "bam"]
NUM_PAIRS   = [16] #[4, 8, 16]
READ_LENS   = [32] #[64, 128]

AIE_FREQ_HZ = 1_250_000_000   # 1.25 GHz (from scsim_config.json)
BITAP_MAX_K = 5
BAM_MIN_MAPQ = 20

# Paths for real Illumina data (illumina and bam modes)
ILLUMINA_FASTQ  = "test_data/illumina_real/NA12878_50k_R1.fastq.gz"
ILLUMINA_FASTA  = "test_data/illumina_real/chr22.fa"
ILLUMINA_BAM    = "test_data/illumina_real/NA12878_chr22_128bp.bam"

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _run(cmd: str, cwd: Path, log_path: Path, label: str):
    """Run shell command, write stdout+stderr to log_path.  Returns (rc, elapsed_s)."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with open(log_path, "w") as lf:
        proc = subprocess.run(
            cmd, shell=True, cwd=str(cwd),
            stdout=lf, stderr=subprocess.STDOUT,
        )
    elapsed = time.time() - t0
    return proc.returncode, elapsed


def _banner(msg: str):
    w = 80
    print(f"\n{'=' * w}")
    print(f"  {msg}")
    print(f"{'=' * w}")


# ─────────────────────────────────────────────────────────────────────────────
# Profile file parsers
# ─────────────────────────────────────────────────────────────────────────────

def parse_profile_funct(path: Path) -> dict:
    """
    Extract from profile_funct_*.txt:
      total_cycles, report_cycles, total_instructions, program_memory_bytes,
      kernel_cycles, kernel_pct, kernel_calls, kernel_name
    """
    result = {
        "total_cycles":        None,
        "report_cycles":       None,
        "total_instructions":  None,
        "program_memory_bytes": None,
        "kernel_cycles":       None,
        "kernel_pct":          None,
        "kernel_calls":        None,
        "kernel_name":         None,
    }
    if not path.exists():
        return result

    text = path.read_text(errors="replace")

    def _int(pattern):
        m = re.search(pattern, text)
        return int(m.group(1).replace(",", "")) if m else None

    result["total_cycles"]         = _int(r"Total cycle count\s*:\s*([\d,]+)")
    result["report_cycles"]        = _int(r"Report cycle count\s*:\s*([\d,]+)")
    result["total_instructions"]   = _int(r"Total instruction count\s*:\s*([\d,]+)")
    result["program_memory_bytes"] = _int(r"Total size in program memory\s*:\s*([\d,]+)")

    # Function table (space-aligned):
    #   calls  cycles_func  %func  min  avg  max  cycles+desc  %  min  avg  max  lowPC  highPC  FunctionName
    # Pick the heaviest non-infrastructure function (= the user kernel)
    skip_funcs = {"main", "_main_init", "memset", "__cxa_finalize", "_fini"}
    func_re = re.compile(
        r"^\s*(\d+)\s+(\d+)\s+([\d.]+)%\s+"   # calls  cycles_func  %func
        r"\d+\s+\d+\s+\d+\s+"                  # min avg max (func)
        r"\d+\s+([\d.]+)%\s+"                  # cycles+desc  %
        r"\d+\s+\d+\s+\d+\s+"                  # min avg max (func+desc)
        r"\d+\s+\d+\s+"                        # low/high PC
        r"(\S+)",                               # first token of function name
        re.MULTILINE,
    )
    best_cycles = -1
    for m in func_re.finditer(text):
        fname = m.group(5)
        if any(fname.startswith(sk) for sk in skip_funcs):
            continue
        cyc = int(m.group(2))
        if cyc > best_cycles:
            best_cycles = cyc
            result["kernel_calls"]  = int(m.group(1))
            result["kernel_cycles"] = cyc
            result["kernel_pct"]    = float(m.group(3))
            result["kernel_name"]   = fname

    return result


# Vector instruction mnemonic prefixes for AIE
_VECTOR_PREFIXES = (
    "VMAC", "VADD", "VSUB", "VMUL", "VMAX", "VMIN",
    "VSRS", "VSHL", "VSHR", "VBROADCAST", "VUPSHIFT",
    "VSHUF", "VABS", "VNEG", "VAND", "VOR", "VXOR",
    "VLDA", "VSTA", "VST", "VLD",
)


def parse_profile_instr(path: Path, kernel_name: Optional[str]) -> dict:
    """
    Extract from profile_instr_*.txt (fields not present in profile_funct):
      report_instr_count, report_instr_coverage,
      kernel_instr_count, kernel_instr_coverage, kernel_prog_mem_bytes,
      total_stall_cycles, nop_count, vector_instr_count, scalar_instr_count
    """
    result = {
        # header fields unique to profile_instr
        "report_instr_count":    None,
        "report_instr_coverage": None,   # % of addresses profiled
        # per-kernel fields from function summary / detail
        "kernel_instr_count":    None,
        "kernel_instr_coverage": None,   # % of kernel code actually executed
        "kernel_prog_mem_bytes": None,
        # instruction-level analysis from kernel detail section
        "total_stall_cycles":    0,
        "nop_count":             0,
        "vector_instr_count":    0,
        "scalar_instr_count":    0,
    }
    if not path.exists():
        return result

    text = path.read_text(errors="replace")

    # ── header ──────────────────────────────────────────────────────────────
    m = re.search(r"Report instruction count\s*:\s*([\d,]+)", text)
    if m:
        result["report_instr_count"] = int(m.group(1).replace(",", ""))

    m = re.search(r"Report instruction coverage\s*:\s*([\d.]+)%", text)
    if m:
        result["report_instr_coverage"] = float(m.group(1))

    # ── function summary table ────────────────────────────────────────────
    # Format:  cycles  % of total  instructions  % of total  % Coverage  FunctionName...
    skip_funcs = {"main", "_main_init", "memset", "__cxa_finalize", "_fini"}
    fsumm_re = re.compile(
        r"^\s+(\d+)\s+([\d.]+)%\s+(\d+)\s+([\d.]+)%\s+([\d.]+)%\s+(\S+)",
        re.MULTILINE,
    )
    best_cyc = -1
    for m in fsumm_re.finditer(text):
        fname = m.group(6)
        if any(fname.startswith(sk) for sk in skip_funcs):
            continue
        cyc = int(m.group(1))
        if cyc > best_cyc:
            best_cyc = cyc
            result["kernel_instr_count"]    = int(m.group(3))
            result["kernel_instr_coverage"] = float(m.group(5))

    # ── kernel function detail section ───────────────────────────────────
    # Find the right detail block by matching kernel_name; fall back to full text
    section = text
    if kernel_name:
        idx = text.find(f"Function detail: {kernel_name}")
        if idx >= 0:
            next_idx = text.find("Function detail:", idx + 1)
            section = text[idx: next_idx if next_idx > 0 else len(text)]

            # Extract size in program memory for this function
            pm = re.search(r"Size in program memory\s*:\s*([\d,]+)", section)
            if pm:
                result["kernel_prog_mem_bytes"] = int(pm.group(1).replace(",", ""))

    # Each instruction line inside the detail block:
    #   <spaces>PC  hex_bytes  assembly  exe_count  cycles  user_count  wait_states  ...
    instr_re = re.compile(
        r"^\s{4,}\d+\s+"              # PC (at least 4 leading spaces to avoid header rows)
        r"(?:[0-9a-f]{2}\s+)+"        # hex bytes
        r"(.+?)\s+"                   # assembly  (non-greedy)
        r"(\d+)\s+"                   # exe_count
        r"\d+\s+"                     # cycles
        r"\d+\s+"                     # user_count
        r"(\d+)",                     # wait_states
        re.MULTILINE,
    )

    for m in instr_re.finditer(section):
        asm_raw = m.group(1).strip()
        exe     = int(m.group(2))
        wait    = int(m.group(3))

        result["total_stall_cycles"] += exe * wait

        # Classify each micro-op in the VLIW bundle (ops separated by ";")
        ops = [op.strip() for op in asm_raw.split(";") if op.strip()]
        for op in ops:
            mnemonic = op.split()[0].upper() if op.split() else ""
            if mnemonic == "NOP":
                result["nop_count"] += exe
            elif any(mnemonic.startswith(vp) for vp in _VECTOR_PREFIXES):
                result["vector_instr_count"] += exe
            else:
                result["scalar_instr_count"] += exe

    return result


def parse_aiesim_log(path: Path) -> dict:
    """Extract PASS/FAIL, sim_wall_time_ns, and matched_count from AIESimulator.log."""
    result = {
        "pass_fail":      "UNKNOWN",
        "sim_time_ns":    None,
        "matched_count":  None,
        "pass_fail_msg":  "",
    }
    if not path.exists():
        return result

    text = path.read_text(errors="replace")

    # PASS: Bitap - 5 scores matched  /  FAIL: ...
    pf_m = re.search(r"^(PASS|FAIL):(.+)$", text, re.MULTILINE)
    if pf_m:
        result["pass_fail"]     = pf_m.group(1)
        result["pass_fail_msg"] = pf_m.group(2).strip()
        cnt = re.search(r"(\d+)\s+scores?\s+matched", pf_m.group(2))
        if cnt:
            result["matched_count"] = int(cnt.group(1))

    # Total Simulation time 92240 ns
    t_m = re.search(r"Total Simulation time\s+([\d.]+)\s+ns", text)
    if t_m:
        result["sim_time_ns"] = float(t_m.group(1))

    return result


def parse_map_report(path: Path) -> dict:
    """Extract memory bank assignments from Map_Report.csv."""
    import csv as _csv
    result = {"mem_banks": "", "num_tiles": 0}
    if not path.exists():
        return result

    banks = []
    tiles = set()
    section = None   # "BUFFER", "PLIO", "CLUSTER"

    with open(path, newline="", errors="replace") as fh:
        for row in _csv.reader(fh):
            if not row:
                continue
            first = row[0].strip()
            if first == "BUFFER":
                section = "BUFFER"
                continue
            if first == "PLIO":
                section = "PLIO"
                continue
            if first == "CLUSTER":
                section = "CLUSTER"
                continue

            if section == "BUFFER" and len(row) >= 6:
                # name, from, to, size, memory_group, bank, offset, dma
                bank      = row[5].strip()
                mem_group = row[4].strip()
                if bank.isdigit():
                    banks.append(f"{mem_group}@bank{bank}")

            elif section == "CLUSTER" and len(row) >= 3:
                # cluster_id, kernels, tile "(col, row)"
                tile = row[2].strip()
                if re.match(r"\(\d+,\s*\d+\)", tile):
                    tiles.add(tile)

    result["mem_banks"] = "; ".join(banks) if banks else "N/A"
    result["num_tiles"] = len(tiles)
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Per-run orchestration
# ─────────────────────────────────────────────────────────────────────────────

def run_one(algo_name, data_algo, kernel_dir, data_type, num_pairs, read_len,
            log_dir: Path, dry_run: bool) -> dict:
    """Run one sweep combination and return a dict of all metrics."""

    aie_dir    = ROOT / "kernels" / kernel_dir / "aie"
    build_dir  = aie_dir / "build.hw"
    sim_out    = build_dir / "aiesimulator_output"
    data_dir   = aie_dir / "data"

    tag = f"{data_type}_{algo_name}_np{num_pairs}_rl{read_len}"
    run_log_dir = log_dir / tag

    row = {
        "data_type":     data_type,
        "algo_name":     algo_name,
        "kernel_dir":    kernel_dir,
        "num_pairs":     num_pairs,
        "read_len":      read_len,
        # will fill in below
        "datagen_status": "",
        "build_status":   "",
        "sim_status":     "",
        "pass_fail":      "N/A",
        "pass_fail_msg":  "",
        "matched_count":  None,
        # ── from profile_funct_24_0.txt ──────────────────────────────────
        "total_cycles":         None,
        "report_cycles":        None,
        "total_instructions":   None,
        "program_memory_bytes": None,   # Total size in program memory
        "kernel_name":          None,
        "kernel_cycles":        None,
        "kernel_pct":           None,   # % of total cycles
        "kernel_calls":         None,
        # ── from profile_instr_24_0.txt ──────────────────────────────────
        "report_instr_count":    None,  # unique to profile_instr header
        "report_instr_coverage": None,  # % instruction addresses profiled
        "kernel_instr_count":    None,  # instructions executed in kernel func
        "kernel_instr_coverage": None,  # % of kernel code paths executed
        "kernel_prog_mem_bytes": None,  # kernel function size in prog memory
        # instruction-level breakdown (from kernel detail section)
        "total_stall_cycles":   None,
        "stall_pct":            None,
        "nop_count":            None,
        "vector_instr_count":   None,
        "scalar_instr_count":   None,
        "issue_rate":           None,   # kernel_instr_count / (kernel_cycles × 7)
        "vec_util_pct":         None,   # vector_instr_count / kernel_instr_count × 100
        # timing
        "sim_time_ns":        None,
        "kernel_time_ns":     None,
        "kernel_time_us":     None,
        # perf
        "mcups":              None,
        "cycles_per_cell":    None,
        # memory
        "mem_banks":          "N/A",
        "num_tiles":          0,
        # notes
        "notes": "",
    }

    label = f"{data_type:10s} | {algo_name:38s} | np={num_pairs:2d} | rl={read_len:3d}"
    print(f"\n  [{label}]")

    if dry_run:
        row["notes"] = "dry-run"
        return row

    # ── 1. Data generation ───────────────────────────────────────────────────
    bitap_flag = f" --bitap-max-k {BITAP_MAX_K}" if data_algo == "bitap" else ""

    if data_type == "bam" and not (ROOT / ILLUMINA_BAM).exists():
        row["datagen_status"] = "SKIP (BAM not found)"
        row["notes"] = f"BAM file missing: {ILLUMINA_BAM}"
        print(f"SKIP (BAM not found)")
        return row

    if data_algo == "swg_interseq":
        # Batched 2-bit interseq format — generated separately from the plio-dp
        # format the other kernels use.
        interseq_common = (
            f" --seq-len {read_len}"
            f" --batch-size {num_pairs}"
            f" --num-batches 2"
            f" --out-dir kernels/{kernel_dir}/aie/data"
            f" --seed 42"
        )
        if data_type == "synthetic":
            datagen_cmd = f"python3 scripts/gen_interseq_data.py{interseq_common}"
        elif data_type == "bam":  # real BAM — extract read/ref windows at true mapping positions
            datagen_cmd = (
                f"python3 scripts/gen_interseq_data.py"
                f" --bam {ILLUMINA_BAM}"
                f" --ref-fasta {ILLUMINA_FASTA}"
                f" --min-mapq {BAM_MIN_MAPQ}"
                f"{interseq_common}"
            )
        else:  # illumina — use real FASTQ + reference FASTA
            datagen_cmd = (
                f"python3 scripts/gen_interseq_data.py"
                f" --fastq {ILLUMINA_FASTQ}"
                f" --ref-fasta {ILLUMINA_FASTA}"
                f"{interseq_common}"
            )
    elif data_type == "synthetic":
        datagen_cmd = (
            f"python3 scripts/generate_test_data.py"
            f" --algorithms {data_algo}"
            f" --out-dir kernels/{kernel_dir}/aie/data"
            f" --num-pairs {num_pairs}"
            f" --ref-len {read_len}"
            f" --read-len {read_len}"
            f"{bitap_flag}"
        )
    elif data_type == "bam":  # real BAM — extract read/ref windows at true mapping positions
        datagen_cmd = (
            f"python3 scripts/generate_test_data.py"
            f" --bam {ILLUMINA_BAM}"
            f" --ref-fasta {ILLUMINA_FASTA}"
            f" --min-mapq {BAM_MIN_MAPQ}"
            f" --algorithms {data_algo}"
            f" --out-dir kernels/{kernel_dir}/aie/data"
            f" --num-pairs {num_pairs}"
            f" --read-len {read_len}"
            f"{bitap_flag}"
        )
    else:  # illumina — use real FASTQ + reference FASTA
        datagen_cmd = (
            f"python3 scripts/generate_test_data.py"
            f" --fastq {ILLUMINA_FASTQ}"
            f" --ref-fasta {ILLUMINA_FASTA}"
            f" --algorithms {data_algo}"
            f" --out-dir kernels/{kernel_dir}/aie/data"
            f" --num-pairs {num_pairs}"
            f" --read-len {read_len}"
            f"{bitap_flag}"
        )

    print(f"    [1/4] Generating {data_type} data ...", end=" ", flush=True)
    rc, elapsed = _run(
        datagen_cmd, ROOT,
        run_log_dir / "datagen.log",
        label="datagen",
    )
    if rc != 0:
        row["datagen_status"] = f"FAIL (rc={rc})"
        row["notes"] = "Data generation failed"
        print(f"FAIL (rc={rc})")
        return row
    row["datagen_status"] = f"OK ({elapsed:.1f}s)"
    print(f"OK ({elapsed:.1f}s)")

    # ── 2. Build ─────────────────────────────────────────────────────────────
    # Use chmod + rm -rf directly instead of make clean to avoid:
    #   "rm: cannot remove 'build.hw/aiesimulator_output': Directory not empty"
    # which occurs when aiesimulator sets read-only permissions on its output dirs.
    print(f"    [2/4] clean + make all ...", end=" ", flush=True)
    clean_and_build = (
        "chmod -R u+rwx build.hw 2>/dev/null; "
        "rm -rf build.hw 2>/dev/null; "
        "make all"
    )
    rc, elapsed = _run(
        clean_and_build,
        aie_dir,
        run_log_dir / "build.log",
        label="build",
    )
    if rc != 0:
        row["build_status"] = f"FAIL (rc={rc}, {elapsed:.0f}s)"
        row["notes"] = "Build failed"
        print(f"FAIL (rc={rc})")
        return row
    row["build_status"] = f"OK ({elapsed:.0f}s)"
    print(f"OK ({elapsed:.0f}s)")

    # ── 3. Simulate ──────────────────────────────────────────────────────────
    print(f"    [3/4] make sim ...", end=" ", flush=True)
    rc, elapsed = _run(
        "make sim",
        aie_dir,
        run_log_dir / "sim.log",
        label="sim",
    )
    if rc != 0:
        row["sim_status"] = f"FAIL (rc={rc}, {elapsed:.0f}s)"
        row["notes"] = "Simulation failed"
        print(f"FAIL (rc={rc})")
        # still try to parse partial results below

    else:
        row["sim_status"] = f"OK ({elapsed:.0f}s)"
        print(f"OK ({elapsed:.0f}s)")

    # ── 4. Parse results ─────────────────────────────────────────────────────
    print(f"    [4/4] Parsing profile outputs ...", end=" ", flush=True)

    # PASS/FAIL from AIESimulator.log
    sim_log = build_dir / "AIESimulator.log"
    log_data = parse_aiesim_log(sim_log)
    row.update({k: log_data[k] for k in log_data})

    # Profile function
    funct_files = sorted(sim_out.glob("profile_funct_*.txt"))
    funct_data  = {}
    if funct_files:
        funct_data = parse_profile_funct(funct_files[0])  # single tile
        row.update({k: funct_data[k] for k in funct_data})

    # Profile instruction
    instr_files = sorted(sim_out.glob("profile_instr_*.txt"))
    if instr_files:
        instr_data = parse_profile_instr(instr_files[0], funct_data.get("kernel_name"))
        row.update({k: instr_data[k] for k in instr_data})

    # Memory map
    map_csv = build_dir / "Map_Report.csv"
    map_data = parse_map_report(map_csv)
    row.update({k: map_data[k] for k in map_data})

    # ── 5. Derived metrics ───────────────────────────────────────────────────
    kernel_cycles = row.get("kernel_cycles")
    total_cycles  = row.get("total_cycles")
    report_cycles = row.get("report_cycles")
    total_instr   = row.get("total_instructions")
    stall_cycles  = row.get("total_stall_cycles", 0) or 0

    if kernel_cycles:
        row["kernel_time_ns"] = kernel_cycles / (AIE_FREQ_HZ / 1e9)   # cycles / (GHz)
        row["kernel_time_us"] = row["kernel_time_ns"] / 1000.0

    VLIW_SLOTS = 7
    kernel_instr = row.get("kernel_instr_count")
    kernel_cyc   = row.get("kernel_cycles")
    if kernel_instr and kernel_cyc and kernel_cyc > 0:
        row["issue_rate"] = round(kernel_instr / (kernel_cyc * VLIW_SLOTS), 4)

    vec_count         = row.get("vector_instr_count") or 0
    kernel_instr_count = row.get("kernel_instr_count")
    if kernel_instr_count and kernel_instr_count > 0:
        # % of kernel instructions that are vector ops
        row["vec_util_pct"] = round(vec_count / kernel_instr_count * 100, 4)

    if report_cycles and stall_cycles is not None:
        row["stall_pct"] = round(stall_cycles / report_cycles * 100, 2) if report_cycles > 0 else None

    # Performance metrics: cells = num_pairs × read_len × read_len
    # (for DP-based algorithms each pair needs read_len × read_len cells)
    num_cells = num_pairs * read_len * read_len
    if kernel_cycles and kernel_cycles > 0:
        row["cycles_per_cell"] = round(kernel_cycles / num_cells, 4)

    if row.get("kernel_time_us") and row["kernel_time_us"] > 0:
        row["mcups"] = round(num_cells / row["kernel_time_us"], 3)

    status_str = row.get("pass_fail", "?")
    matched    = row.get("matched_count")
    mcups_str  = f"{row['mcups']:.1f} MCUPS" if row.get("mcups") else "N/A"
    print(f"{status_str}  |  {mcups_str}")

    return row


# ─────────────────────────────────────────────────────────────────────────────
# CSV output
# ─────────────────────────────────────────────────────────────────────────────

# (header, row-dict key)  — order defines CSV column order
COLUMNS = [
    ("data_type",            "data_type"),
    ("algo_name",            "algo_name"),
    ("kernel_dir",           "kernel_dir"),
    ("num_pairs",            "num_pairs"),
    ("read_len",             "read_len"),
    ("pass_fail",            "pass_fail"),
    ("matched_count",        "matched_count"),
    ("datagen_status",       "datagen_status"),
    ("build_status",         "build_status"),
    ("sim_status",           "sim_status"),
    # ── profile_funct_24_0.txt ───────────────────────────────────────────────
    ("total_cycles",         "total_cycles"),
    ("report_cycles",        "report_cycles"),
    ("total_instructions",   "total_instructions"),
    ("program_memory_bytes", "program_memory_bytes"),
    # ── kernel function breakdown ─────────────────────────────────────────────
    ("kernel_name",          "kernel_name"),
    ("kernel_cycles",        "kernel_cycles"),
    ("kernel_pct",           "kernel_pct"),
    ("kernel_calls",         "kernel_calls"),
    # ── profile_instr_24_0.txt ───────────────────────────────────────────────
    ("report_instr_count",   "report_instr_count"),
    ("report_instr_coverage","report_instr_coverage"),
    ("kernel_instr_count",   "kernel_instr_count"),
    ("kernel_instr_coverage","kernel_instr_coverage"),
    ("kernel_prog_mem_bytes","kernel_prog_mem_bytes"),
    ("vliw_issue_rate",      "issue_rate"),
    ("stall_cycles",         "total_stall_cycles"),
    ("stall_pct",            "stall_pct"),
    ("nop_count",            "nop_count"),
    ("vector_instr_count",   "vector_instr_count"),
    ("scalar_instr_count",   "scalar_instr_count"),
    ("vec_util_pct",         "vec_util_pct"),
    # ── timing / performance ─────────────────────────────────────────────────
    ("sim_wall_time_ns",     "sim_time_ns"),
    ("kernel_time_ns",       "kernel_time_ns"),
    ("kernel_time_us",       "kernel_time_us"),
    ("mcups",                "mcups"),
    ("cycles_per_cell",      "cycles_per_cell"),
    # ── memory ───────────────────────────────────────────────────────────────
    ("mem_banks",            "mem_banks"),
    ("num_tiles",            "num_tiles"),
    ("notes",                "notes"),
]

CSV_HEADERS = [h for h, _ in COLUMNS]


def write_csv(rows: list, output_path: Path):
    """Write all rows to CSV, overwriting the file each time (safe incremental saves)."""
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(CSV_HEADERS)
        for row in rows:
            writer.writerow([row.get(key, "") for _, key in COLUMNS])
    print(f"\n  CSV saved → {output_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="AIESIM sweep for genomics kernels")
    parser.add_argument("--dry-run",  action="store_true",
                        help="Print combinations without running anything")
    parser.add_argument("--output",   default="results_aiesim_sweep.csv",
                        help="Output CSV filename (default: results_aiesim_sweep.csv)")
    parser.add_argument("--algo",     nargs="+", choices=list(ALGO_KERNEL_MAP.keys()),
                        help="Restrict to specific algo names")
    parser.add_argument("--data-type", nargs="+", choices=DATA_TYPES,
                        help="Restrict to specific data types")
    parser.add_argument("--num-pairs", nargs="+", type=int,
                        help="Override num_pairs list (default: 8 16 32)")
    parser.add_argument("--read-len",  nargs="+", type=int,
                        help="Override read_len list  (default: 64 128)")
    args = parser.parse_args()

    algos      = {k: v for k, v in ALGO_KERNEL_MAP.items()
                  if args.algo is None or k in args.algo}
    dtypes     = args.data_type or DATA_TYPES
    num_pairs  = args.num_pairs or NUM_PAIRS
    read_lens  = args.read_len  or READ_LENS

    total_runs = len(dtypes) * len(algos) * len(num_pairs) * len(read_lens)
    output_path = ROOT / args.output

    log_dir = ROOT / "logs_aiesim_sweep" / datetime.now().strftime("%Y%m%d_%H%M%S")

    # ── BAM pre-flight check ─────────────────────────────────────────────────
    if "bam" in dtypes:
        bam_abs  = ROOT / ILLUMINA_BAM
        bai_abs  = Path(str(bam_abs) + ".bai")
        fasta_abs = ROOT / ILLUMINA_FASTA
        missing = []
        if not bam_abs.exists():
            missing.append(f"BAM   : {bam_abs}")
        if not bai_abs.exists():
            missing.append(f"BAI   : {bai_abs}")
        if not fasta_abs.exists():
            missing.append(f"FASTA : {fasta_abs}")
        if missing:
            print("\nERROR: BAM mode requires the following files (not found):")
            for m in missing:
                print(f"  {m}")
            print("\nGenerate them with:")
            print(f"  bwa index {ILLUMINA_FASTA}")
            print(f"  bwa mem {ILLUMINA_FASTA} {ILLUMINA_FASTQ} \\")
            print(f"      | samtools sort -o {ILLUMINA_BAM}")
            print(f"  samtools index {ILLUMINA_BAM}")
            raise SystemExit(1)

    _banner(f"AIESIM Sweep — {total_runs} combinations")
    print(f"  Data types : {dtypes}")
    print(f"  Algos      : {list(algos.keys())}")
    print(f"  Num pairs  : {num_pairs}")
    print(f"  Read lens  : {read_lens}")
    print(f"  Output     : {output_path}")
    print(f"  Logs       : {log_dir}")
    if "bam" in dtypes:
        print(f"  BAM        : {ILLUMINA_BAM}  (min_mapq={BAM_MIN_MAPQ})")
    if args.dry_run:
        print("  [DRY RUN — no builds will be executed]")

    all_rows = []
    run_n    = 0

    for dt in dtypes:
        _banner(f"Data type: {dt}")
        for algo_name, (data_algo, kernel_dir) in algos.items():
            for np in num_pairs:
                for rl in read_lens:
                    run_n += 1
                    print(f"\n  Run {run_n}/{total_runs}", end="")
                    row = run_one(
                        algo_name, data_algo, kernel_dir,
                        dt, np, rl,
                        log_dir, args.dry_run,
                    )
                    all_rows.append(row)
                    # Save incrementally after each run
                    write_csv(all_rows, output_path)

    _banner("Sweep complete")
    print(f"  Total runs    : {total_runs}")
    passed = sum(1 for r in all_rows if r.get("pass_fail") == "PASS")
    failed = sum(1 for r in all_rows if r.get("pass_fail") == "FAIL")
    print(f"  PASS          : {passed}")
    print(f"  FAIL / other  : {total_runs - passed}")
    print(f"  Results       : {output_path}")
    print(f"  Logs          : {log_dir}")


if __name__ == "__main__":
    main()
