#!/usr/bin/env python3
"""
gen_pl_utilization_report.py -- Record PL/AIE resource utilization (%) for one
genomics_server algo build.

Resource utilization is a property of the hardware build (algo/kernel +
ROWS/COLS/BUF_SIZE), not of the runtime dataset (synthetic/illumina/bam) --
feeding different data doesn't change LUT/FF/BRAM/URAM/DSP usage. So this
script produces ONE row per algo (not per data_type), reading whichever
Vivado post-place utilization report is currently on disk.

kernels/genomics_server/_x/reports/link/imp/impl_1_full_util_placed.rpt is
overwritten by every `make link`/`make package` run, so it only ever reflects
the algo that was most recently built. Run this script once right after each
algo's hardware build to accumulate a full table across algos:

  make -C kernels/genomics_server aie pl link package KERNEL=wfa ROWS=8 COLS=50 BUF_SIZE=2048
  python3 scripts/gen_pl_utilization_report.py --algo genomics_server_wfa

  make -C kernels/genomics_server aie pl link package KERNEL=bitap ROWS=8 COLS=50 BUF_SIZE=2048
  python3 scripts/gen_pl_utilization_report.py --algo genomics_server_bitap_k5
  ...

Repeated runs for the same --algo overwrite that algo's row (idempotent).

Usage:
  python3 scripts/gen_pl_utilization_report.py --algo genomics_server_wfa
  python3 scripts/gen_pl_utilization_report.py --list   # show recorded rows
"""

import argparse
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REPORT = ROOT / "kernels/genomics_server/_x/reports/link/imp/impl_1_full_util_placed.rpt"
KERNEL_CONFIG_H = ROOT / "kernels/genomics_server/aie/build.hw/kernel_config.h"
CSV_OUT = ROOT / "kernels/genomics_server/results_pl_utilization.csv"

# Fixed config this report is meaningful for (row=8 col=50 buf_size=2048).
# num_iters / refs_per_batch are datagen/runtime params -- they don't affect
# the hardware build, so they aren't checked here.
EXPECTED_ROWS = 1
EXPECTED_COLS = 1
EXPECTED_BUF_SIZE = 2048

# KERNEL_DEF (from kernel_config.h) -> algo name, for the unambiguous cases.
# bitap is ambiguous (k5 vs k100 are different hardware builds via MAX_K in
# bitap_kernel.cc / bitap_kernel_v1.cc) so it's intentionally left out --
# --algo must be given explicitly for bitap builds.
KERNEL_DEF_HINT = {
    "USE_WFA":           "genomics_server_wfa",
    "USE_STRIPED_SWG":   "genomics_server_striped_swg",
    "USE_INTERSEQ_SWG":  "genomics_server_interseq",
    "USE_SCALAR_SWG":    "genomics_server_scalar_swg",
}

METRICS = [
    ("lut_pct",   "CLB LUTs"),
    ("ff_pct",    "Registers"),
    ("bram_pct",  "Block RAM Tile"),
    ("uram_pct",  "URAM"),
    ("dsp_pct",   "DSP Slices"),
]

CSV_FIELDS = ["algo", "rows", "cols", "buf_size", "lut_pct", "ff_pct", "bram_pct", "uram_pct", "dsp_pct"]


def parse_kernel_config(path: Path) -> dict:
    if not path.exists():
        return {}
    text = path.read_text()
    cfg = {}
    for key in ("NUM_ROWS", "NUM_COLS", "BUF_SIZE", "NUM_ITERS"):
        m = re.search(rf"#define {key}\s+(\d+)", text)
        if m:
            cfg[key] = int(m.group(1))
    m = re.search(r"#define (USE_\w+)", text)
    if m:
        cfg["KERNEL_DEF"] = m.group(1)
    return cfg


def parse_utilization(path: Path) -> dict:
    if not path.exists():
        print(f"ERROR: report not found: {path}", file=sys.stderr)
        print("  Run `make link` / `make package` for the target algo first.", file=sys.stderr)
        sys.exit(1)

    lines = path.read_text(errors="replace").splitlines()
    result = {}
    for field, label in METRICS:
        found = None
        for line in lines:
            if "|" not in line:
                continue
            parts = [p.strip() for p in line.strip().split("|")]
            if len(parts) >= 7 and parts[1] == label:
                found = float(parts[6])
                break
        if found is None:
            print(f"WARNING: could not find row '{label}' in {path}", file=sys.stderr)
        result[field] = found
    return result


def load_csv(path: Path) -> dict:
    rows = {}
    if path.exists():
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                rows[row["algo"]] = row
    return rows


def write_csv(path: Path, rows: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for algo in sorted(rows):
            writer.writerow(rows[algo])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--algo", help="Algo name to record this report under (e.g. genomics_server_wfa)")
    parser.add_argument("--report", default=str(REPORT), help="Path to impl_1_full_util_placed.rpt")
    parser.add_argument("--force", action="store_true",
                         help="Record even if kernel_config.h doesn't match the expected ROWS/COLS/BUF_SIZE")
    parser.add_argument("--list", action="store_true", help="Print the current CSV and exit")
    args = parser.parse_args()

    rows = load_csv(CSV_OUT)

    if args.list:
        if not rows:
            print(f"No rows recorded yet ({CSV_OUT} doesn't exist).")
            return
        w = max(len(a) for a in rows) + 2
        print(f"{'algo':<{w}}{'LUT%':>8}{'FF%':>8}{'BRAM%':>8}{'URAM%':>8}{'DSP%':>8}")
        for algo, row in sorted(rows.items()):
            print(f"{algo:<{w}}{row['lut_pct']:>8}{row['ff_pct']:>8}{row['bram_pct']:>8}"
                  f"{row['uram_pct']:>8}{row['dsp_pct']:>8}")
        print(f"\n({CSV_OUT})")
        return

    if not args.algo:
        cfg = parse_kernel_config(KERNEL_CONFIG_H)
        hint = KERNEL_DEF_HINT.get(cfg.get("KERNEL_DEF"))
        if hint:
            print(f"--algo not given; auto-detected '{hint}' from {KERNEL_CONFIG_H}")
            args.algo = hint
        else:
            print("ERROR: --algo is required (auto-detect failed or is ambiguous, e.g. bitap "
                  "k5 vs k100 both show USE_BITAP).", file=sys.stderr)
            sys.exit(1)

    cfg = parse_kernel_config(KERNEL_CONFIG_H)
    mismatch = []
    if cfg.get("NUM_ROWS") != EXPECTED_ROWS:
        mismatch.append(f"NUM_ROWS={cfg.get('NUM_ROWS')} (expected {EXPECTED_ROWS})")
    if cfg.get("NUM_COLS") != EXPECTED_COLS:
        mismatch.append(f"NUM_COLS={cfg.get('NUM_COLS')} (expected {EXPECTED_COLS})")
    if cfg.get("BUF_SIZE") != EXPECTED_BUF_SIZE:
        mismatch.append(f"BUF_SIZE={cfg.get('BUF_SIZE')} (expected {EXPECTED_BUF_SIZE})")
    if mismatch and not args.force:
        print(f"ERROR: {KERNEL_CONFIG_H} config doesn't match the expected build:\n  "
              + "\n  ".join(mismatch), file=sys.stderr)
        print("  Rebuild with ROWS=8 COLS=50 BUF_SIZE=2048, or pass --force to record anyway.",
              file=sys.stderr)
        sys.exit(1)

    util = parse_utilization(Path(args.report))

    row = {
        "algo": args.algo,
        "rows": cfg.get("NUM_ROWS", "N/A"),
        "cols": cfg.get("NUM_COLS", "N/A"),
        "buf_size": cfg.get("BUF_SIZE", "N/A"),
        **util,
    }
    rows[args.algo] = row
    write_csv(CSV_OUT, rows)

    print(f"Recorded {args.algo}: "
          f"LUT={util['lut_pct']}%  FF={util['ff_pct']}%  BRAM={util['bram_pct']}%  "
          f"URAM={util['uram_pct']}%  DSP={util['dsp_pct']}%")
    print(f"-> {CSV_OUT}")


if __name__ == "__main__":
    main()
