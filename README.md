<img src="https://github.com/Juanitovh/CACTUS/assets/53839626/1440936f-f3aa-453f-bed6-8afc28b26b2c" width="200" height="200" />
<img src="https://github.com/Juanitovh/CACTUS/assets/53839626/9a43c5be-980e-4e98-a8dd-4dfe7d05e192" width="200" height="200" />
<img src="https://github.com/Juanitovh/CACTUS/assets/53839626/ccbc15b8-ec17-45af-8b63-25331edf0582" width="200" height="200" />


# CACTUS
CACTUS: A Computational Framework for Generating Realistic White Matter Microstructure Substrates

CACTUS webpage [http://cactus.epfl.ch/](http://cactus.epfl.ch/) 

Frontiers [Paper](https://www.frontiersin.org/articles/10.3389/fninf.2023.1208073/full)



# CACTUS Installation Guide in bash Ubunt

Follow these steps to install and set up **CACTUS** on Ubuntu.

## Step 1: Update the System
Update your package manager and upgrade installed packages:
```bash
sudo apt update && sudo apt upgrade -y
```
## Step 2: Install C++ Compiler, OpenMP, python3, pip, and python-venv
```bash 
sudo apt install -y build-essential g++ libomp-dev
sudo apt install -y python3 python3-pip python3-venv
```

## Step 3 Create a Virtual Environment and source it
```bash
python3 -m venv cactus_env
source cactus_env/bin/activate
```

## Step 3 Clone the repo
```bash
git clone https://github.com/Juanitovh/CACTUS.git
cd CACTUS
```

## Step 4:  Install dependencies and compile cpp code
```bash
pip install -r requirements_cactus.txt
bash compile_cpp_code.sh
```

## Step 5:Optional (recomended) Add CACTUS/cactus_scripts to your path

```bash
echo "export PATH=\"$(pwd)/cactus_scripts:\$PATH\"" >> ~/.bashrc
source ~/.bashrc
python update_paths.py

```

## Step 6: Optional, software to visualize the meshes, or your favorite mesh viewer
```bash
sudo apt install meshlab
```

# Ready to run a quick example


## CACTUS Usage Overview

The CACTUS pipeline is now organized into separate wrapper scripts, each handling a specific step of the process outlined in the paper. This structure makes it easy to run the key processes in sequence. 

### Wrapper Scripts:
1. **`1_wrapper_initialization.py`**: Initializes the substrate using the configuration provided in `cactus_single.txt`.
2. **`2_wrapper_optimization.py`**: Performs joint fibre optimization.
3. **`3_wrapper_FRG_mesh.py`**: Executes the Fibre Radial Growth (FRG) algorithm.

### Configuration File:
The file **`cactus_single.txt`** contains all the necessary morphological parameters and hyperparameters to create the substrates, ensuring each step is properly configured and ready for processing.



## Tutorial, toy examples. Can run in personal computer
To begin the tutorial begin by copying "cactus_single.txt" to a folder of your choice. 


## First Step
This step creates 
*.init files, that contains the fibre initializations of the substrate.

```bash
1_wrapper_initialization.py -config_file cactus_single.txt
```



### Substep (optionall create a mesh of the init file)
We can use and include script to quickly generate a mesh using cylinder.
Then we can use meshlab to open it
```bash
quick_mesh_cactus.py tutorial_single_00000.init
meshlab test2.ply
```




## Second step Step: Global joint optimization

This step will start optimize the substrate using gradient descent. It creates auxiliary files to save different time points of the optimized substrate. 
The *.partial files are version of the substrate save every 50 iterations.

If the optimization procedure converges, It creates a optimized_final.txt file containing the final version of the substrate
```bash
 2_wrapper_optimization.py -config_file cactus_single.txt
```

### Substep (optionall create a mesh of the partial file or the final file)
We can use and include script to quickly generate a mesh using cylinder.
Then we can use meshlab to open it
For any of the partials.
```bash
quick_mesh_cactus.py tutorial_single_00001/optimized_00100.partial
meshlab test2.ply
```

For the final optimized version
```bash
quick_mesh_cactus.py tutorial_single_00001/optimized_final.txt
meshlab test2.ply
```


### Substep (optional) 
You can check the status of the optimization by running the following command. 
It allows for real time updates to visualize the status the optimization process.
```bash
optim_loger.py -folder tutorial_single_00000/
```

## Third Step: Fibre Radial Growth
This step is the last step of the pipeline. It will create a mesh of the optimized substrate and apply the FRG algorithm to it.
It needs the optimized_final.txt file to run. Also, it's quite computationally expensive.



### FRG: Select a Case and Substep

To use the CACTUS pipeline, you need to choose a **case** to run and specify the appropriate **substep**. 

- The **case** determines  which strands to process.
- The **substep** defines whether you are growing or meshing the fibres.

### 1. **Choosing a Case**:
- **test**: Run a small batch of fibres to check hyperparameters and ensure correct processing.
- **missing**: Process fibres that are missing or were wrongly computed.
- **all**: Process 100% of the fibres (not recommended on a single laptop due to high computational demand).

### 2. **Choosing a Substep**:
- **growth**: Grow fibres in discrete space (slow process, needed to prepare the structure).
- **mesh**: Mesh the fibres (fast process, done after growth).

### Example Usage:

1. **Run the "growth" substep in the "test" case**:
```bash
python 3_wrapper_FRG_mesh.py -config_file cactus_single.txt -substep growth -run_case test -file tutorial_single_00000.init
```



2. Run the "mesh" substep in the "test" case:

This command will grow the fibres for a small test batch using the cactus_single.txt configuration file and the tutorial_single_00000.init file for input.

```bash
python 3_wrapper_FRG_mesh.py -config_file cactus_single.txt -substep mesh -run_case test -file tutorial_single_00000.init
```



3. Optional: Try the "missing" case and check the results:
```bash
python 3_wrapper_FRG_mesh.py -config_file cactus_single.txt -substep growth -run_case missing -file tutorial_single_00000.init
```


### Output

After running the FRG algorithm, the output will be saved in the following directories:

- **`meshes/`**: Contains all the outputs
- **`meshes/pickles/`**: Stores FRG metadata in the form of sparse compressed matrices. **It's heavy on memory**
- **`/meshes/simulations/`**: Contains meshes prepared for simulations.

### Visualizing the Meshes

To visualize the generated meshes using **MeshLab**, you can use the following command:

```bash
meshlab tutorial_single_00000/meshes/simulations/*.ply
```





# Done
***

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
