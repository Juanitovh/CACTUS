#!/usr/bin/env python3

import argparse
import os
import numpy as np

# Import custom libraries
from libraries import RotatingStrand as RS
from libraries import WattsonFunction as WF
from libraries import ReadWriteStrands as RWS
from libraries import ReadConfigurations as RC
from libraries.CactusPaths import cactus_paths
from libraries.CactusMath import exponential_linspace_int



from libraries import ascii_art 
ascii_art.display_title()

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Second Step Algorithm")
    parser.add_argument(
        "-config_file", 
        type=str, 
        help="Configuration file to initialize several voxel settings. Overrides all parameters.", 
        default=None
    )
    parser.add_argument(
        "-file", 
        type=str, 
        help="*.init file to be processed instead of whole directory.",
        default=None
    )
    parser.add_argument(
    "-run_case",
    help="""Select a CASE to run:
    - test: Run a small batch of fibres to check hyperparameters and ensure proper processing.
    - missing: Grow/mesh all fibres that are missing or were wrongly computed.
    - all: Grow/mesh 100% of the fibres. NOT RECOMMENDED on a single laptop.""",
    choices=["test", "missing", "all"],
    default="test"
    )
    parser.add_argument(
        '-substep', help="""Select one of the two substeps:
        growth: Growing in discrete space, it's slow.
        mesh: Meshing the fibres, it's fast.""",
        choices=["growth", "mesh"]
    )
    return parser.parse_args()


def main():
    """Main function for the script."""
    ascii_art.print_message("Third Step Algorithm ")

    # Parse arguments
    args = parse_arguments()

    # Read configuration file
    if args.config_file ==None:
        print("Error! give me a cactus config_file")
        exit()
    else :
        cactus_args = RC.read_config_file(args.config_file)


    if args.file is not None:
        if not args.file.endswith(".init"):
            print(f"Error!!! the file: {args.file} should be a *.init file")
            print("Exiting...")
            exit()

        files = [args.file]
        print("Processing file: " , files)
    else:
        files = RC.get_files_with_extension(".init")
        print(f"Discovered files: {files}")



    current_experiment_folder = os.getcwd()
    for i, file_name in enumerate(files):
        # Skip experiments as specified
        n=  RC.count_number_streamlines(file_name)

        folder= file_name.split(".")[0]
        os.chdir(current_experiment_folder)
        os.chdir(folder)




        #create file to run axons 
        missing_axon_file =  "0_missing_axon_file.txt" 
        #final_path_missing = os.path.join(folder , missing_axon_file ) 

        if os.path.isfile(missing_axon_file):
            os.remove(missing_axon_file)


        ### case 1: grow some axons for testing

        command_check = None
        if args.run_case == "test":
            missing_strands = exponential_linspace_int(0,n-1 , 5)
            print("")
            print("-"*50)
            print("Strands to be processed for testing: ")
            print(missing_strands)
            print("-"*50)
            print("")
            np.savetxt(missing_axon_file, missing_strands , fmt="%d")

        elif args.run_case == "missing":
            if args.substep == "growth":
                script_check_pickles= os.path.join(cactus_paths.cactus_code , "check_pickles_whole.py" )
                command_check =f"python {script_check_pickles} -file {file_name} -outfile {missing_axon_file} "
            elif args.substep == "mesh":
                script_check_meshes= os.path.join(cactus_paths.cactus_code , "check_meshes_whole.py" )
                command_check =f"python {script_check_meshes} -file {file_name} -outfile {missing_axon_file} "
            else:
                print("Error!!! Select a substep")
                exit()
        elif args.run_case == "all":
            if os.path.isfile(missing_axon_file):
                os.remove(missing_axon_file)
            missing_axon_file = "all"

        if command_check != None:
            ascii_art.print_message("Checking missing axons")
            print("-"*50)
            os.system
            os.system(command_check)
        

        ascii_art.print_message(" Finished checking which strands to run \n the next step next step in the processing")


        if args.substep == "growth":
            script_run_pickles= os.path.join(cactus_paths.cactus_code , "meta_grid.py" )
            command_run =f"python {script_run_pickles} -file {file_name} -missing_axon_file  {missing_axon_file}  -iterations {cactus_args.growth_iterations} -grid_size {cactus_args.grid_size}"
        elif args.substep == "mesh":
            script_run_meshes= os.path.join(cactus_paths.cactus_code , "bake_mesh_pickle.py" )
            command_run =f"python {script_run_meshes} -file {file_name} -missing_axon_file  {missing_axon_file} -colorless {cactus_args.colorless} -inn_out {cactus_args.inn_out} -sim_vol {cactus_args.sim_vol} -n_erode {cactus_args.n_erode} -g_ratio {cactus_args.g_ratio} \n"
        else:
            ascii_art.print_error_message("Error!!! Select a run_case")
            exit()


        print("#" * 20)
        print(command_run)
        print("#" * 20)
        os.system(command_run)



        print("#" * 20)
        ascii_art.print_message("Done!")

if __name__ == "__main__":
    main()
