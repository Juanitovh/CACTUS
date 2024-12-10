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

## Step 5: Add CACTUS/cactus_scripts to your path
```bash
echo "export PATH=\"$(pwd)/cactus_scripts:\$PATH\"" >> ~/.bashrc
source ~/.bashrc
python update_paths.py
python 
```

## Step 6: Optional, software to visualize the meshes
```bash
sudo apt install meshlab
```

# Ready to run a quick example







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
