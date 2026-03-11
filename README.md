<div align="center">

# CACTUS

### Computational Axonal Configurator for Tailored and Ultradense Substrates

[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/downloads/)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-00599C.svg)](https://en.cppreference.com/w/cpp/20)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Paper](https://img.shields.io/badge/Frontiers-10.3389%2Ffninf.2023.1208073-red.svg)](https://www.frontiersin.org/articles/10.3389/fninf.2023.1208073/full)

*A computational framework for generating realistic white matter microstructure substrates for Monte-Carlo diffusion simulations and DW-MRI model validation.*

**[Website](http://cactus.epfl.ch/)** &bull; **[Paper](https://www.frontiersin.org/articles/10.3389/fninf.2023.1208073/full)** &bull; **[Tutorial Data](https://drive.google.com/drive/folders/1G6rz6WjFr7Z5Ii7P16ymfE9KHYm_YVHL?usp=sharing)**

<br>

<table>
  <tr>
    <td><img src="https://github.com/Juanitovh/CACTUS/assets/53839626/1440936f-f3aa-453f-bed6-8afc28b26b2c" width="200" height="200" /></td>
    <td><img src="https://cactus.epfl.ch/images/icvf2.gif" width="200" height="200" /></td>
    <td><img src="https://cactus.epfl.ch/images/ezgif.com-resize.gif" width="200" height="200" /></td>
  </tr>
  <tr>
    <td><img src="https://cactus.epfl.ch/images/cross_002.png" width="200" height="200" /></td>
    <td><img src="https://cactus.epfl.ch/images/ezgif.com-resize1.gif" width="200" height="200" /></td>
    <td><img src="https://cactus.epfl.ch/images/high_pack.png" width="200" height="200" /></td>
  </tr>
  <tr>
    <td><img src="https://cactus.epfl.ch/images/animation.gif" width="200" height="200" /></td>
    <td><img src="https://cactus.epfl.ch/images/ezgif.com-resize.gif" width="200" height="200" /></td>
    <td><img src="https://cactus.epfl.ch/images/disco_mesh.png" width="200" height="200" /></td>
  </tr>
</table>

</div>

<br>

## Overview

CACTUS generates realistic synthetic white matter substrates with high packing densities, large voxel sizes, and diverse fibre configurations. The framework implements a three-stage pipeline:

1. **Initialization** &mdash; Place fibres on a grid using Watson distribution for angular dispersion and gamma distribution for radii
2. **Optimization** &mdash; Refine fibre arrangement via Adagrad gradient descent (C++ with OpenMP) to minimize overlaps and reach target packing density
3. **Growth & Meshing** &mdash; Apply the Fibre Radial Growth (FRG) algorithm with marching cubes to produce 3D PLY meshes

The generated substrates are designed for Monte-Carlo diffusion simulators such as [MCDC](https://github.com/jonhrafe/MCDC_Simulator_public).

---

## Branches

| Branch | Description |
|--------|-------------|
| `main` | Original version from the published paper |
| `master` | **v1.5** &mdash; Significant refactoring, pip-installable package with unified CLI |

---

## Installation

### Prerequisites

- **Python** >= 3.9
- **g++** with C++20 and OpenMP support

```bash
# Ubuntu / Debian
sudo apt install build-essential g++ libomp-dev
```

### Quick start

```bash
git clone https://github.com/Juanitovh/CACTUS.git
cd CACTUS
make          # Compiles C++ optimizer + installs Python package
```

This makes the `cactus1-substrates` CLI available system-wide.

<details>
<summary><b>Manual install (advanced)</b></summary>

```bash
make build-cpp       # Compile C++ optimizer only
pip install -e .     # Install Python package in editable mode
```

</details>

---

## Usage

CACTUS provides a unified CLI with subcommands for each pipeline stage.

### Stage 1 &mdash; Initialize fibre placements

Creates `.init` files with initial fibre positions based on morphological parameters.

```bash
cactus1-substrates init -config_file config_files/cactus_single.txt
```

### Stage 2 &mdash; Global joint optimization

Refines fibre arrangement using gradient descent. Produces `optimized_final.txt`.

```bash
cactus1-substrates optimize -config_file config_files/cactus_single.txt
```

### Stage 3 &mdash; Fibre radial growth & meshing

Applies the FRG algorithm and generates PLY meshes via marching cubes.

```bash
# Growth substep (computationally intensive)
cactus1-substrates grow -config_file config_files/cactus_single.txt -substep growth -run_case test

# Mesh substep (fast, follows growth)
cactus1-substrates grow -config_file config_files/cactus_single.txt -substep mesh -run_case test
```

> **Run case options:** `test` (small batch) &bull; `missing` (only incomplete fibres) &bull; `all` (full dataset)

### Quick mesh generation

Generate a mesh directly from a strand file:

```bash
cactus1-substrates quick-mesh -file optimized_final.txt
cactus1-substrates quick-mesh -file optimized_final.txt --parallel -output output.ply
```

### Utilities

```bash
# Monitor optimization progress in real time
cactus1-substrates monitor -folder my_experiment/

# Convert NFG format to CACTUS format
cactus1-substrates convert -folder nfg_strands/ -outfile output.txt
```

---

## Configuration

Configuration files use space-separated key-value pairs. Multiple values on a single line create parameter sweeps.

| Config file | Use case |
|-------------|----------|
| `config_files/cactus_single.txt` | Single fibre bundle substrates |
| `config_files/cactus_crossing.txt` | Crossing fibre bundle substrates |

---

## Output

| Extension | Stage | Description |
|-----------|:-----:|-------------|
| `.init` | 1 | Initial fibre placements |
| `.partial` | 2 | Optimization checkpoints (every 50 iterations) |
| `optimized_final.txt` | 2 | Final optimized fibre configuration |
| `.bz2` pickles | 3 | Compressed voxel-space metadata |
| `.ply` | 3 | 3D triangle meshes |

### Visualizing meshes

```bash
meshlab tutorial_single_00000/meshes/simulations/*.ply
```

---

## Architecture

```
src/cactus1_substrate/
├── cli.py               # Unified CLI entry point
├── paths.py             # Dynamic path resolution
├── core/                # Library modules
│   ├── cactus_math, rotating_strand, wattson_function, ...
├── pipeline/            # Three-stage pipeline wrappers
│   ├── initialization, optimization, growth_mesh
├── workers/             # Parallelized heavy computation
│   ├── meta_grid, bake_mesh_pickle, grid_initialization, ...
└── tools/               # Standalone utilities
    ├── quick_mesh, optim_logger, nfg2cactus, ...
```

---

## Data & Resources

<details>
<summary><b>24 Hours of DIFFUSION Around the World (ISMRM Tutorial)</b></summary>

<br>

**[Download tutorial data](https://drive.google.com/drive/folders/1G6rz6WjFr7Z5Ii7P16ymfE9KHYm_YVHL?usp=sharing)**

5 substrates with mean fibre radii of 0.20, 0.25, 0.35, 0.60, and 0.75 &mu;m:
- Voxel size: (40 &mu;m)&sup3;
- Fibre radii range: 0.15 &ndash; 2.00 &mu;m
- ICVF: 90%

<div align="center">
<img src="https://github.com/Juanitovh/CACTUS/assets/53839626/42853e86-b037-4837-a00a-5cff18f5f684" width="130" />
<img src="https://github.com/Juanitovh/CACTUS/assets/53839626/ad92a0ae-2f47-4dff-8c45-c6918e6cae44" width="130" />
<img src="https://github.com/Juanitovh/CACTUS/assets/53839626/e14dfff2-9d81-412d-9280-51f5cda4f6cf" width="130" />
<img src="https://github.com/Juanitovh/CACTUS/assets/53839626/d7997b8c-1901-4baa-8b78-df3e6b9493f2" width="130" />
<img src="https://github.com/Juanitovh/CACTUS/assets/53839626/d2d49e28-bac0-465b-9189-eac73540f441" width="130" />
</div>

</details>

<details>
<summary><b>Paper data &mdash; Synthetic twins & DW-MRI simulations</b></summary>

<br>

**[Download paper data](https://drive.google.com/drive/folders/1S2cdEin0uO91FJpUNGTH_I7ZdVKcNClu?usp=sharing)**

<div align="center">
<img src="https://github.com/Juanitovh/CACTUS/assets/53839626/643d1d58-2b3e-4ebb-badf-21c3d01655e0" width="60%" />
</div>

</details>

---

## Citation

If you use CACTUS in your research, please cite:

```bibtex
@article{villarreal2023cactus,
  title     = {CACTUS: A Computational Framework for Generating Realistic
               White Matter Microstructure Substrates},
  author    = {Villarreal-Haro, Juan Luis and others},
  journal   = {Frontiers in Neuroinformatics},
  year      = {2023},
  doi       = {10.3389/fninf.2023.1208073}
}
```

---

<div align="center">
<sub>Developed at <b>EPFL</b> &bull; Signal Processing Laboratory (LTS5)</sub>
</div>
