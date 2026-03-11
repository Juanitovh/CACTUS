# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CACTUS (Computational Axonal Configurator for Tailored and Ultradense Substrates) generates realistic synthetic white matter microstructure substrates for Monte-Carlo diffusion simulations and DW-MRI model validation. It is a hybrid Python/C++ scientific computing framework, packaged as an installable Python package with a CLI.

## Setup Commands

```bash
# Full build: compile C++ optimizer + install Python package
make

# Or separately:
make build-cpp           # Compile C++ optimizer
pip install -e .         # Install Python package in editable mode
```

## Running the Pipeline

The framework provides a single CLI entry point `cactus1-substrates` with subcommands:

```bash
# Stage 1: Initialize fibre placements → produces .init files
cactus1-substrates init -config_file config_files/cactus_single.txt

# Stage 2: Global joint optimization (runs C++ optimizer) → produces optimized_final.txt
cactus1-substrates optimize -config_file config_files/cactus_single.txt

# Stage 3: Fibre radial growth & meshing → produces .ply mesh files
cactus1-substrates grow -config_file config_files/cactus_single.txt -substep growth -run_case test
cactus1-substrates grow -config_file config_files/cactus_single.txt -substep mesh -run_case test

# Quick mesh from strand file
cactus1-substrates quick-mesh -file my_strands.txt
cactus1-substrates quick-mesh -file my_strands.txt --parallel -output output.ply

# Monitor optimization
cactus1-substrates monitor -folder my_experiment/

# Convert NFG → CACTUS format
cactus1-substrates convert -folder nfg_strands/ -outfile output.txt
```

Use `cactus_single.txt` for single-bundle substrates and `cactus_crossing.txt` for crossing fibre bundles.

## Architecture

### Package structure (`src/cactus1_substrate/`)

```
src/cactus1_substrate/
├── __init__.py          # Package version
├── __main__.py          # python -m cactus1_substrate
├── cli.py               # argparse CLI with subcommands
├── paths.py             # Dynamic path resolution (auto-finds optimizer binary)
├── core/                # Library modules (snake_case)
│   ├── ascii_art.py, cactus_math.py, read_configurations.py,
│   ├── read_write_strands.py, rotating_strand.py, slurm.py,
│   ├── wattson_function.py, wavy_radii.py
├── pipeline/            # Three-stage pipeline wrappers
│   ├── initialization.py, optimization.py, growth_mesh.py
├── workers/             # Subprocess-invoked heavy scripts
│   ├── grid_initialization.py, meta_grid.py, bake_mesh_pickle.py,
│   ├── check_pickles.py, check_meshes.py
└── tools/               # Standalone utility commands
    ├── quick_mesh.py, quick_mesh_parallel.py, optim_logger.py,
    ├── nfg2cactus.py, paste_mesh_dataset.py
```

### Pipeline stages

1. **Initialization** (`pipeline/initialization.py` → `workers/grid_initialization.py`): Places fibres on a grid using Watson distribution for dispersion, gamma distribution for radii. Outputs `.init` files.

2. **Optimization** (`pipeline/optimization.py` → `joint_fibre_optimizer`): C++ executable using Adagrad gradient descent with OpenMP parallelization. Minimizes fibre overlaps to achieve target packing density. Saves `.partial` checkpoints every 50 iterations.

3. **Growth & Meshing** (`pipeline/growth_mesh.py` → `workers/meta_grid.py`, `tools/quick_mesh.py`): Voxel-space discretization, radial growth, and marching cubes to produce PLY meshes. Supports parallel batch processing.

### C++ optimizer (`source_fibre_optimization/`)

- Built with `make build-cpp` (uses Makefile)
- Compiled with: `g++ -O3 -Ofast -std=c++20 -fopenmp`
- Output binary: `source_fibre_optimization/joint_fibre_optimizer`

### Path resolution (`paths.py`)

The optimizer binary is found automatically via fallback chain:
1. `CACTUS_OPTIMIZER_PATH` environment variable
2. `shutil.which("joint_fibre_optimizer")` (binary on PATH)
3. `<repo_root>/source_fibre_optimization/joint_fibre_optimizer` (Makefile build location)

No manual path configuration needed (replaces old `update_paths.py`).

### Key output artifacts

| Extension | Stage | Description |
|-----------|-------|-------------|
| `.init` | 1 | Initial fibre placements |
| `.partial` | 2 | Optimization checkpoints |
| `optimized_final.txt` | 2 | Final optimized fibre configuration |
| `.ply` | 3 | 3D mesh files |
| `.bz2` pickles | 3 | Compressed mesh metadata |

## Important Notes

- No automated test suite exists. Validation is done via intermediate file inspection and MeshLab visualization.
- Numba JIT is used extensively in numerical routines — first invocations will be slow due to compilation.
- The config file format uses space-separated key-value pairs; multiple values on one line create parameter sweeps.
- Workers are invoked as `python -m cactus1_substrate.workers.<module>` from the pipeline wrappers.
