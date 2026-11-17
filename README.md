# VersalAlign

Sequence-alignment kernels (WFA, Bitap, striped SWG, inter-sequence SWG, scalar SWG) on the AMD Versal AI Engine array, targeting the **VCK5000** data-center card.

This file covers **how to compile** the designs. To regenerate the numbers reported in the paper, see [reproduce_results.md](reproduce_results.md).

## Repository layout

```text
genomics-alignment/
├── kernels/
│   ├── common/            # shared headers (dna_utils.hh, parameters.hh)
│   ├── genomics_server/   # main design: R x C AIE array + resident PL servers (used in the paper)
│   │   ├── aie/           #   AIE graph; selects the compute kernel via KERNEL=...
│   │   ├── pl/            #   HLS kernels dist_server / collect_server (+ gen_servers.py)
│   │   ├── hw_link/       #   v++ link configs (8x50: config.cfg, 1x1: config_1.cfg)
│   │   ├── host/          #   XRT host application
│   │   └── scripts/       #   gen_multi_data.py (input + golden data for the array)
│   ├── wfa/  bitap/  swg_interseq/  swg_scalar/   # single-tile kernels: aie/ pl/ hw_link/ host/
│   └── sw_striped/        # single-tile kernel: aie/ only (simulation)
├── scripts/               # data generation, sweep, power and utilization scripts
└── test_data/illumina_real/
```

The `genomics_server` AIE graph compiles the per-algorithm source files from `kernels/<algo>/aie/src/` directly, so those directories are needed even if you only build `genomics_server`.

## 1. Prerequisites

| Item | Version used |
|---|---|
| AMD Vitis (includes `aiecompiler`, `v++`, `aiesimulator`, `x86simulator`) | **2023.1** |
| XRT | version that matches the VCK5000 platform/shell |
| Platform | `xilinx_vck5000_gen4x8_qdma_2_202220_1` |
| Host compiler | `g++` with C++17 |
| Python | 3.x; `pysam` for BAM input |

Set up the environment in every shell before building:

```bash
source <VITIS_INSTALL>/Vitis/2023.1/settings64.sh
source /opt/xilinx/xrt/setup.sh

# Check the tools and platform
which aiecompiler v++
platforminfo -p xilinx_vck5000_gen4x8_qdma_2_202220_1 --json="file"   # must print a .xpfm path
echo $XILINX_VITIS $XILINX_XRT                                         # used by the Makefiles
```

If `platforminfo` doesn't find the platform, add its install directory to `PLATFORM_REPO_PATHS`.

## 2. Generate input data (needed for simulation and board runs)

The AIE graphs read their PLIO inputs and golden outputs from `aie/data/`. That directory isn't checked in, so generate it first. Run these commands from `genomics-alignment/`:

```bash
cd genomics-alignment

# genomics_server array, synthetic data:
#   <out_dir> <rows> <cols> <refs_per_row> <num_iters> <buf_size> [seed]
python3 kernels/genomics_server/scripts/gen_multi_data.py \
    kernels/genomics_server/aie/data 8 50 200 1 2048 42 --kernel wfa --read-len 128
#   --kernel is wfa | sw | bitap  (use sw for striped_swg, interseq_swg and scalar_swg)
#   add  --bam <file.bam> --ref-fasta <ref.fa>  to use real reads instead

# Single-tile kernels:
python3 scripts/generate_test_data.py --algorithms wfa --out-dir kernels/wfa/aie/data \
    --num-pairs 16 --ref-len 32 --read-len 32
python3 scripts/gen_interseq_data.py --bam test_data/illumina_real/NA12878_chr22_128bp.bam \
    --ref-fasta test_data/illumina_real/chr22.fa --seq-len 32 --num-batches 2 \
    --out-dir kernels/swg_interseq/aie/data --min-mapq 20 --seed 42 --batch-size 16
```

[reproduce_results.md](reproduce_results.md) §1 explains how to download the real NA12878 / chr22 dataset and build the BAM file.

## 3. Build the main design: `genomics_server`

### 3.1 Choose the configuration

The array size **must be the same in four places**. If they differ, the link step fails, or the design hangs or fails on the board.

| What | Where | Default |
|---|---|---|
| AIE graph size and kernel | `kernels/genomics_server/aie/Makefile`: `KERNEL`, `ROWS`, `COLS`, `BUF_SIZE`, `SCORE_BUF`, `ITERS` (or pass them on the `make` command line) | `wfa`, 8, 50, 2048, 64, 1 |
| PL server kernels | `kernels/genomics_server/pl/gen_servers.py`: `NUM_ROWS`, `NUM_COLS`, `BUF_SIZE`, `SCORE_BUF`, `RING` | 8, 50, 2048, 64, 2 |
| Stream connectivity | `kernels/genomics_server/hw_link/config.cfg` (one `stream_connect` per row input and per column output) | 8 inputs / 50 outputs |
| Host | `#define`s at the top of `kernels/genomics_server/host/host_genomics_server.cpp`: `NUM_ROWS`, `NUM_COLS`, `BUF_SIZE`, `SCORE_BUF`, `NUM_ITERS`, `REF_PER_ROW`, `RING` | 8, 50, 2048, 64, 1000, 200, 2 |

Valid values for `KERNEL`:

| `KERNEL=` | Algorithm | AIE stack size |
|---|---|---|
| `wfa` | Wavefront alignment (gap-affine edit distance) | 8000 |
| `bitap` | Bitap / Wu–Manber, `MAX_K` set in `kernels/bitap/aie/src/bitap_kernel.cc` (default 5) | 6000 |
| `striped_swg` | Striped Smith–Waterman–Gotoh | 6000 |
| `interseq_swg` | Inter-sequence SWG | 6000 |
| `scalar_swg` | Scalar SWG | 6000 |
| anything else | dummy kernel (for testing data movement only) | 4000 |

To change the array size (for example to 1x1):

```bash
cd genomics-alignment/kernels/genomics_server

# 1) PL: edit NUM_ROWS / NUM_COLS in pl/gen_servers.py, then regenerate the server sources
(cd pl && python3 gen_servers.py)          # rewrites pl/dist_server.cpp and pl/collect_server.cpp

# 2) Link config: the top-level Makefile always uses hw_link/config.cfg
cp hw_link/config.cfg hw_link/config_8x50.cfg.bak
cp hw_link/config_1.cfg hw_link/config.cfg  # 1x1 (8x50 is the default config.cfg / config_400.cfg)
# Other sizes need a config with in_0..in_{R-1} and out_0..out_{C-1}.

# 3) Host: edit NUM_ROWS / NUM_COLS / REF_PER_ROW in host/host_genomics_server.cpp

# 4) AIE: pass ROWS / COLS when building (next section)
```

### 3.2 Compile

Run from `genomics-alignment/kernels/genomics_server/`. Variables given on the command line are passed down to the `aie/` sub-make.

```bash
cd genomics-alignment/kernels/genomics_server

make clean                                    # required whenever KERNEL/ROWS/COLS change (see note)
make all KERNEL=wfa ROWS=8 COLS=50            # aie + pl + link + package + host
```

Or run the steps one at a time:

| Target | What it runs | Output |
|---|---|---|
| `make aie` | `aiecompiler --target=hw` on `aie/src/array_graph_v3m.cpp` | `aie/build.hw/libadf.a` |
| `make pl` | `v++ -c` (HLS, 300 MHz) for `dist_server` and `collect_server` | `pl/build.hw/*.xo` |
| `make link` | `v++ -l` with `hw_link/config.cfg` | `build.hw/genomics_server.xsa` |
| `make package` | `v++ -p` | `build.hw/genomics_server.xclbin` |
| `make host` | `g++ -std=c++17` + XRT | `host/genomics_server` |

> **Note:** `aie/build.hw/kernel_config.h` (which holds `NUM_ROWS`, `NUM_COLS`, `USE_<KERNEL>` and so on) is generated **only if it doesn't already exist**. If you change `KERNEL`, `ROWS`, `COLS` or `ITERS` without running `make clean` (or `make -C aie clean`), the previous configuration is built again.

The hardware link and package steps take several hours for the 8x50 array.

### 3.3 Name the xclbin per algorithm

`scripts/run_sweep.sh` looks for the xclbins in `genomics-alignment/` under these names:

```bash
# after building with KERNEL=wfa
cp build.hw/genomics_server.xclbin ../../genomics_server_wfa.xclbin
```

| `KERNEL=` | Expected xclbin name |
|---|---|
| `wfa` | `genomics_server_wfa.xclbin` |
| `bitap` (`MAX_K=5`) | `genomics_server_bitap_k5.xclbin` |
| `striped_swg` | `genomics_server_striped_swg.xclbin` |
| `interseq_swg` | `genomics_server_interseq.xclbin` |
| `scalar_swg` | `genomics_server_scalar_swg.xclbin` |

Repeat sections 3.2–3.3 once for each algorithm.

### 3.4 Simulate before building for hardware (optional)

Simulation is much faster than a hardware build. Generate `aie/data/` first (section 2) with the same rows and columns.

```bash
make x86sim KERNEL=wfa ROWS=1 COLS=1    # functional simulation (x86simulator)
make aiesim KERNEL=wfa ROWS=1 COLS=1    # cycle-approximate simulation (aiesimulator --profile)
```

The graph's `main()` compares `aiesimulator_output/output_col_<c>.txt` (or `x86simulator_output/...`) against `aie/data/golden_<c>.txt` and prints PASS or FAIL. `make clean` before switching between `x86sim` and `hw`, because both use `kernel_config.h`.

### 3.5 Run on the board

```bash
cd genomics-alignment/kernels/genomics_server/host
./genomics_server <path/to/xclbin> ../aie/data [bdf]   # default bdf: 0000:82:00.1
```

Use `xbutil examine` to find your card's BDF. The run prints `PASS`/`FAIL` and the throughput in GCUPS.

## 4. Build the single-tile kernels (optional)

`wfa`, `bitap`, `swg_interseq` and `swg_scalar` each come with a 1-tile design: `aie/`, `pl/` (`pl_mm2s` / `pl_s2mm` data movers), `hw_link/config.cfg` and `host/`. `sw_striped` has only `aie/`. These directories have **no top-level Makefile**, so you run the link and package steps yourself:

```bash
cd genomics-alignment/kernels/wfa          # or bitap, swg_interseq, swg_scalar
PLATFORM=xilinx_vck5000_gen4x8_qdma_2_202220_1

make -C aie                                # -> aie/build.hw/libadf.a   (TARGET=hw by default)
make -C pl                                 # -> pl/build.hw/pl_mm2s.xo, pl_s2mm.xo
mkdir -p build.hw
v++ -l -t hw --platform $PLATFORM --config hw_link/config.cfg \
    aie/build.hw/libadf.a pl/build.hw/pl_mm2s.xo pl/build.hw/pl_s2mm.xo \
    -o build.hw/wfa.xsa
v++ -p -t hw --platform $PLATFORM build.hw/wfa.xsa aie/build.hw/libadf.a \
    -o build.hw/wfa.xclbin
make -C host                               # -> host/wfa.exe (bitap.exe, ...)

./host/wfa.exe build.hw/wfa.xclbin aie/data [debug]
```

AIE simulation only (this is what `scripts/aiesim_sweep.py` runs for each algorithm):

```bash
make -C aie                  # TARGET=hw  -> aiesimulator
make -C aie sim
# or, for the functional simulator:
make -C aie clean && make -C aie TARGET=x86sim && make -C aie sim TARGET=x86sim
```

## 5. Clean

```bash
make -C genomics-alignment/kernels/genomics_server clean   # aie + pl + host + build.hw
make -C genomics-alignment/kernels/<algo>/aie clean        # single-tile kernels: clean each subdir
```

Build outputs (`build.*`, `*.xo`, `*.xsa`, `*.xclbin`, logs, simulator output) are listed in `genomics-alignment/.gitignore`.
