# CACTUS

## A Computational Framework for Generating Realistic White Matter Microstructure Substrates

<div align="center">
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

CACTUS (Computational Axonal Configurator for Tailored and Ultradense Substrates) is a computational workflow for generating realistic synthetic white matter substrates. It enables the creation of complex tissue structures with high packing densities, large voxel sizes, and diverse fibre configurations. The generated substrates are ideal for Monte-Carlo diffusion simulations, aiding in the validation of diffusion-weighted MRI (DW-MRI) models.

For more details, refer to:
- **[CACTUS Website](http://cactus.epfl.ch/)**
- **[Published Paper](https://www.frontiersin.org/articles/10.3389/fninf.2023.1208073/full)**

> **Branches:** The `main` branch contains the original version of CACTUS as published in the paper. The `master` branch is version 1.5 of the CACTUS library, featuring significant changes: proper Python packaging (`pip install`), a unified CLI (`cactus1-substrates`), Makefile-based C++ builds, snake_case module naming, and dynamic path resolution.

---

## Installation

### Prerequisites
- Python >= 3.10
- g++ with C++20 and OpenMP support (`sudo apt install build-essential g++ libomp-dev`)

### Install

```bash
git clone https://github.com/Juanitovh/CACTUS.git
cd CACTUS
make
```

This compiles the C++ optimizer and installs the Python package in editable mode. The `cactus1-substrates` CLI command becomes available.

### Alternative: manual install

```bash
# Compile C++ optimizer
make build-cpp

# Install Python package
pip install -e .
```

---

## Usage

CACTUS provides a single CLI entry point with subcommands for each pipeline stage:

### Step 1: Initialize fibre placements

Creates `.init` files with initial fibre positions based on morphological parameters.

```bash
cactus1-substrates init -config_file config_files/cactus_single.txt
```

### Step 2: Global joint optimization

Refines fibre arrangement using gradient descent. Produces `optimized_final.txt`.

```bash
cactus1-substrates optimize -config_file config_files/cactus_single.txt
```

### Step 3: Fibre radial growth and meshing

Applies the FRG algorithm and generates PLY meshes.

```bash
# Growth substep (slow, computationally intensive)
cactus1-substrates grow -config_file config_files/cactus_single.txt -substep growth -run_case test

# Mesh substep (fast, follows growth)
cactus1-substrates grow -config_file config_files/cactus_single.txt -substep mesh -run_case test
```

Run case options: `test` (small batch), `missing` (only incomplete fibres), `all` (full dataset).

### Quick mesh generation

Generate a mesh directly from a strand file:

```bash
cactus1-substrates quick-mesh -file my_strands.txt
cactus1-substrates quick-mesh -file my_strands.txt --parallel -output output.ply
```

### Monitor optimization progress

```bash
cactus1-substrates monitor -folder my_experiment/
```

### Convert NFG format to CACTUS format

```bash
cactus1-substrates convert -folder nfg_strands/ -outfile output.txt
```

---

## Configuration

Use `config_files/cactus_single.txt` for single-bundle substrates and `config_files/cactus_crossing.txt` for crossing fibre bundles. The config file format uses space-separated key-value pairs; multiple values on one line create parameter sweeps.

---

## Output Directories
- **`meshes/`**: Stores all generated meshes
- **`meshes/pickles/`**: Contains compressed metadata
- **`meshes/simulations/`**: Holds simulation-ready meshes

### Visualizing the Meshes
```bash
meshlab tutorial_single_00000/meshes/simulations/*.ply
```

---

## References

#### 24 Hours of DIFFUSION Around the World: ISMRM
**Tutorial data:**
  - Small meshes toy example: [data](https://drive.google.com/drive/folders/1G6rz6WjFr7Z5Ii7P16ymfE9KHYm_YVHL?usp=sharing)

Details:
- 5 substrates: mean fibre radii of 0.20 um, 0.25 um, 0.35 um , 0.60 um , 0.75 um
- Voxel size (40 um)^3
- Minimum fibre radii 0.15 um
- Maximum fibre radii 2.00 um
- ICVF 90%
- Config file to Run simulations with the [MCDC](https://github.com/jonhrafe/MCDC_Simulator_public)


<img src="https://github.com/Juanitovh/CACTUS/assets/53839626/42853e86-b037-4837-a00a-5cff18f5f684" width="150" height="150" />
<img src="https://github.com/Juanitovh/CACTUS/assets/53839626/ad92a0ae-2f47-4dff-8c45-c6918e6cae44" width="150" height="150" />
<img src="https://github.com/Juanitovh/CACTUS/assets/53839626/e14dfff2-9d81-412d-9280-51f5cda4f6cf" width="150" height="150" />
<img src="https://github.com/Juanitovh/CACTUS/assets/53839626/d7997b8c-1901-4baa-8b78-df3e6b9493f2" width="150" height="150" />
<img src="https://github.com/Juanitovh/CACTUS/assets/53839626/d2d49e28-bac0-465b-9189-eac73540f441" width="150" height="150" />


***

#### CACTUS
**Paper data**
  - Synthetic twins meshes and DW-MRI simulations: [data](https://drive.google.com/drive/folders/1S2cdEin0uO91FJpUNGTH_I7ZdVKcNClu?usp=sharing)

Frontiers [Paper](https://www.frontiersin.org/articles/10.3389/fninf.2023.1208073/full)
<img src="https://github.com/Juanitovh/CACTUS/assets/53839626/643d1d58-2b3e-4ebb-badf-21c3d01655e0" width=50% />



***
