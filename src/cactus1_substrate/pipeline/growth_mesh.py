#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys
import numpy as np

# Import custom libraries
from cactus1_substrate.core import rotating_strand as RS
from cactus1_substrate.core import wattson_function as WF
from cactus1_substrate.core import read_write_strands as RWS
from cactus1_substrate.core import read_configurations as RC
from cactus1_substrate.paths import cactus_paths
from cactus1_substrate.core.cactus_math import exponential_linspace_int



from cactus1_substrate.core import ascii_art

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
        '-substep', help="Select one of the two substeps: 1) growth: Growing in discrete space, it's slow. 2) mesh: Meshing the fibres, it's fast.",
        choices=["growth", "mesh"],
        type=str
    )
    parser.add_argument(
    "-run_case",
    help="CASE to run. 1) test: Run a small batch of fibres to check hyperparameters and ensure proper processing. 2) missing: process all fibres that are missing or were wrongly computed. 3) all: process all the fibres. NOT RECOMMENDED on a single laptop.",
    choices=["test", "missing", "all"],
    type=str,
    )
    return parser.parse_args()


def main(args=None):
    """Main function for the script."""
    ascii_art.display_title()
    ascii_art.print_message("Third Step Algorithm ")

    # Parse arguments
    if args is None:
        args = parse_arguments()

    # Read configuration file
    if args.config_file is None:
        print("Error! give me a cactus config_file")
        exit()
    else:
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
                command_check = [sys.executable, "-m", "cactus1_substrate.workers.check_pickles",
                                 "-file", file_name, "-outfile", missing_axon_file]
            elif args.substep == "mesh":
                command_check = [sys.executable, "-m", "cactus1_substrate.workers.check_meshes",
                                 "-file", file_name, "-outfile", missing_axon_file]
            else:
                print("Error!!! Select a substep")
                exit()
        elif args.run_case == "all":
            if os.path.isfile(missing_axon_file):
                os.remove(missing_axon_file)
            missing_axon_file = "all"

        if command_check is not None:
            ascii_art.print_message("Checking missing axons")
            print("-"*50)
            subprocess.run(command_check, check=True)


        ascii_art.print_message(" Finished checking which strands to run \n the next step next step in the processing")


        if args.substep == "growth":
            command_run = [sys.executable, "-m", "cactus1_substrate.workers.meta_grid",
                           "-file", file_name, "-missing_axon_file", missing_axon_file,
                           "-iterations", str(cactus_args.growth_iterations),
                           "-grid_size", str(cactus_args.grid_size)]
        elif args.substep == "mesh":
            command_run = [sys.executable, "-m", "cactus1_substrate.workers.bake_mesh_pickle",
                           "-file", file_name, "-missing_axon_file", missing_axon_file,
                           "-colorless", str(cactus_args.colorless),
                           "-inn_out", str(cactus_args.inn_out),
                           "-sim_vol", str(cactus_args.sim_vol),
                           "-n_erode", str(cactus_args.n_erode),
                           "-g_ratio", str(cactus_args.g_ratio)]
        else:
            ascii_art.print_error_message("Error!!! Select a run_case")
            exit()


        print("#" * 20)
        print(command_run)
        print("#" * 20)
        subprocess.run(command_run, check=True)



        print("#" * 20)
        ascii_art.print_message("Done!")


if __name__ == "__main__":
    main()
