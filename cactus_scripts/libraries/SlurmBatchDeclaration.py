
import subprocess

dict_job_type = {
    "meta": 12,
    "bake": 3,
    "doom": 50, ## only applies for intra
    }

import os 
class jobSlurm:
    def __init__(self , current_path=None , folder=None ):

        self.path_experiment = current_path
        self.sub_experiment = folder
        self.choosen_job = None

    def set_subexperiment_path(self, sub_experiment=None):
        self.sub_experiment = sub_experiment
        self.final_path = os.path.join(self.path_experiment , self.sub_experiment)


    def sbatchDoom(self ,args):
        file_out = f"doom_signals.sh"
        f = open(file_out , 'w')
        texto = (
                f"#!/bin/bash \n"
                f"#SBATCH --chdir {self.path_experiment}\n"
                f"#SBATCH --nodes 1 \n"
                f"#SBATCH --ntasks 1 \n"
                f"#SBATCH --cpus-per-task 15\n"
                f"#SBATCH --mem  14000 \n"
                f"#SBATCH --time 06:00:00 \n"
                f"#SBATCH --qos free \n"

                f"python3 ~/strand_optimization/Doom_doom.py -pos {args.intra_extra} -folder_experiment {self.sub_experiment} -n_batchs {args.n_batchs} -batch_id {args.batch_id} -folder_signals {args.folder_signals} -icvf {args.icvf} -propagator {args.propagator}  -subgrid_size {args.subgrid_size}   -subdivisions_number {args.subdivisions_number}  -n_erode {args.n_erode} \n  "
        )
        f.write(texto)
        f.close()

        print("#############################################")
        subprocess.run(f"pwd" , shell=True)
        subprocess.run(f"sbatch {file_out}" , shell=True)
        subprocess.run(f"echo SI SE PUDO {file_out}  in {self.final_path}"    , shell=True)
        print("#############################################")
        return file_out

    def sbatchMetaGrid(self, args):
        file_out = f"meta_grid.sh"
        f = open(file_out , 'w')
        texto = (
                f"#!/bin/bash \n"
                f"#SBATCH --chdir {self.final_path}\n"
                f"#SBATCH --nodes 1 \n"
                f"#SBATCH --ntasks 1 \n"
                f"#SBATCH --cpus-per-task 10\n"
                f"#SBATCH --mem 80000 \n"
                f"#SBATCH --time 06:00:00 \n"
                f"#SBATCH --qos free \n"
                #f"python3 /home/jvillarr/strand_optimization/meta_grid.py -grid_size .2 -batch_id {i} -n_batch {n_splits} -file {file}  -iterations 5"
                f"python3 /home/jvillarr/strand_optimization/meta_grid_broken_tree.py -grid_size {args.grid_size} -batch_id {args.batch_id} -n_batch {args.n_batchs} -file final.txt  -iterations {args.growth_iterations} -missing_axons {args.missing_axons}"
        )
        f.write(texto)
        f.close()
        
        print("#############################################")
        subprocess.run(f"sbatch {file_out}" , shell=True)
        subprocess.run(f"echo SI SE PUDO {file_out}  in {self.final_path}"    , shell=True)
        print("#############################################")
        return file_out


    def sbatchBakeMesh(self ,args):

        file_out = f"bake_mesh.sh"
        f = open(file_out , 'w')
        texto = (
            f"#!/bin/bash \n"
                #f"#SBATCH --chdir /scratch/jvillarr/{folder}\n"
                f"#SBATCH --chdir {self.final_path}\n"
                f"#SBATCH --nodes 1 \n"
                f"#SBATCH --ntasks 1 \n"
                f"#SBATCH --cpus-per-task 5\n"
                f"#SBATCH --mem 10000 \n"
                f"#SBATCH --time 04:00:00 \n"
                f"#SBATCH --qos free \n"
                #f"python3 /home/jvillarr/strand_optimization/bake_mesh_pickle.py -batch_id {i} -n_batch {n_splits} -file {file} -colorless -inn_out outer -sim_vol volume -n_erode 0\n"
                #f"python3 /home/jvillarr/strand_optimization/bake_mesh_pickle.py -batch_id {i} -n_batch {n_splits} -file {file} -colorless -inn_out outer -sim_vol simulations -n_erode 0\n"
                f"python3 /home/jvillarr/strand_optimization/bake_mesh_pickle_broken.py -batch_id {args.batch_id} -n_batch {args.n_batchs} -file final.txt -colorless {args.colorless} -inn_out {args.inn_out} -sim_vol {args.sim_vol} -n_erode {args.n_erode}  -missing_axons {args.missing_axons} \n"
               # f"python3 /home/jvillarr/strand_optimization/bake_mesh_pickle.py       -batch_id {i} -n_batch {n_splits} -file {file} -inn_out inner -sim_vol volume -n_erode 1\n"
        )
        f.write(texto)
        f.close()


        print("#############################################")
        subprocess.run(f"echo SI SE PUDO {file_out}  in {self.final_path}"    , shell=True)
        subprocess.run(f"sbatch {file_out}" , shell=True)
        print("#############################################")
        return file_out






