<img src="https://github.com/Juanitovh/CACTUS/assets/53839626/1440936f-f3aa-453f-bed6-8afc28b26b2c" width="200" height="200" />
<img src="https://github.com/Juanitovh/CACTUS/assets/53839626/9a43c5be-980e-4e98-a8dd-4dfe7d05e192" width="200" height="200" />
<img src="https://github.com/Juanitovh/CACTUS/assets/53839626/ccbc15b8-ec17-45af-8b63-25331edf0582" width="200" height="200" />
<img src="https://cactus.epfl.ch/images/icvf2.gif" width="200" height="200" />
<img src="https://cactus.epfl.ch/images/ezgif.com-resize1.gif" width="200" height="200" />
<img src="https://cactus.epfl.ch/images/high_pack.png" width="200" height="200" />
<img src="https://cactus.epfl.ch/images/animation.gif" width="200" height="200" />
<img src="https://cactus.epfl.ch/images/ezgif.com-resize.gif" width="200" height="200" />
<img src="https://cactus.epfl.ch/images/disco_mesh.png" width="200" height="200" />


# CACTUS

## A Computational Framework for Generating Realistic White Matter Microstructure Substrates

CACTUS (Computational Axonal Configurator for Tailored and Ultradense Substrates) is a computational workflow for generating realistic synthetic white matter substrates. It enables the creation of complex tissue structures with high packing densities, large voxel sizes, and diverse fibre configurations. The generated substrates are ideal for Monte-Carlo diffusion simulations, aiding in the validation of diffusion-weighted MRI (DW-MRI) models.

For more details, refer to:
- **[CACTUS Website](http://cactus.epfl.ch/)**
- **[Published Paper](https://www.frontiersin.org/articles/10.3389/fninf.2023.1208073/full)**

---

# Installation Guide (Ubuntu)

Follow these steps to install and set up **CACTUS** on Ubuntu.

## Step 1: Update System Packages
```bash
sudo apt update && sudo apt upgrade -y
```

## Step 2: Install Dependencies
Install the required C++ compiler, OpenMP, Python, and related libraries:
```bash
sudo apt install -y build-essential g++ libomp-dev
sudo apt install -y python3 python3-pip python3-venv
```

## Step 3: Create and Activate a Virtual Environment
```bash
python3 -m venv cactus_env
source cactus_env/bin/activate
```

## Step 4: Clone the Repository
```bash
git clone https://github.com/Juanitovh/CACTUS.git
cd CACTUS
```

## Step 5: Install Python Dependencies and Compile C++ Code
```bash
pip install -r requirements_cactus.txt
bash compile_cpp_code.sh
```

## Step 6 (Optional): Add CACTUS Scripts to System Path
This step simplifies running CACTUS commands from any location.
```bash
echo "export PATH=\"$(pwd)/cactus_scripts:\$PATH\"" >> ~/.bashrc
source ~/.bashrc
python update_paths.py
```

## Step 7 (Optional): Install a Mesh Viewer (e.g., MeshLab)
```bash
sudo apt install meshlab
```

---

# Running a Quick Example

## Overview of CACTUS Pipeline
CACTUS consists of several wrapper scripts that execute different steps of the pipeline sequentially.

### Wrapper Scripts:
1. **`1_wrapper_initialization.py`**: Initializes the substrate based on parameters in `cactus_single.txt`.
2. **`2_wrapper_optimization.py`**: Optimizes the fibres using gradient descent.
3. **`3_wrapper_FRG_mesh.py`**: Performs Fibre Radial Growth (FRG) and meshing.

### Configuration File:
The **`cactus_single.txt`** file contains all necessary morphological parameters and hyperparameters for substrate generation.

---

# Tutorial: Running CACTUS Step-by-Step

## Step 1: Initialize the Substrate
This step creates `.init` files containing initial fibre placements.
```bash
1_wrapper_initialization.py -config_file cactus_single.txt
```

### (Optional) Generate a Mesh from Initialization File
To visualize the initialized substrate:
```bash
quick_mesh_cactus.py tutorial_single_00000.init
meshlab test2.ply
```

---

## Step 2: Global Joint Optimization
This step refines the fibre arrangement using gradient descent. Intermediate states are saved as `.partial` files every 50 iterations. The final optimized substrate is stored in `optimized_final.txt`.
```bash
2_wrapper_optimization.py -config_file cactus_single.txt
```

### (Optional) Generate a Mesh from Optimization Output
For an intermediate `.partial` file:
```bash
quick_mesh_cactus.py tutorial_single_00001/optimized_00100.partial
meshlab test2.ply
```
For the final optimized substrate:
```bash
quick_mesh_cactus.py tutorial_single_00001/optimized_final.txt
meshlab test2.ply
```

### (Optional) Monitor Optimization Progress in Real Time
```bash
optim_loger.py -folder tutorial_single_00000/
```

---

## Step 3: Fibre Radial Growth (FRG)
This final step applies the FRG algorithm to generate detailed fibre structures. It requires `optimized_final.txt` and is computationally intensive.

### Selecting a Case and Substep
**Run Case Options:**
- **test**: Small batch of fibres (recommended for debugging)
- **missing**: Processes only missing or incorrect fibres
- **all**: Full dataset (requires significant computational power)

**Substep Options:**
- **growth**: Expands fibres in discrete space (slow but necessary)
- **mesh**: Converts fibres into a mesh format (fast, follows growth step)

### Example Usage
Run the **growth** substep on a small test batch:
```bash
3_wrapper_FRG_mesh.py -config_file cactus_single.txt -substep growth -run_case test -file tutorial_single_00000.init
```

Run the **mesh** substep for the test case:
```bash
3_wrapper_FRG_mesh.py -config_file cactus_single.txt -substep mesh -run_case test -file tutorial_single_00000.init
```

Process missing fibres:
```bash
3_wrapper_FRG_mesh.py -config_file cactus_single.txt -substep growth -run_case missing -file tutorial_single_00000.init
```

---

## Output Directories
- **`meshes/`**: Stores all generated meshes
- **`meshes/pickles/`**: Contains compressed metadata (memory-intensive)
- **`meshes/simulations/`**: Holds simulation-ready meshes

### Visualizing the Meshes
To open the generated mesh files in MeshLab:
```bash
meshlab tutorial_single_00000/meshes/simulations/*.ply
```

---

# CACTUS is Ready to Use!
Your CACTUS pipeline is now set up and ready for generating realistic white matter microstructure substrates. 🎉


#### 24 Hours of DIFFUSION Around the World: ISMRM
**Tutorial data:**
  -Small meshes toy example: [data](https://drive.google.com/drive/folders/1G6rz6WjFr7Z5Ii7P16ymfE9KHYm_YVHL?usp=sharing)

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
