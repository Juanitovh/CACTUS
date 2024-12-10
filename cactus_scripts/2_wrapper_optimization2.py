#!/usr/bin/env python3

import argparse
import os
import numpy as np

# Import custom libraries
from libraries import RotatingStrand as RS
from libraries import WattsonFunction as WF
from libraries import ReadWriteStrands as RWS
from libraries import ReadConfigurations as RC

from libraries import ascii_art 

from libraries.CactusPaths import cactus_paths

from libraries import ascii_art 

ascii_art.display_title()
ascii_art.print_message("Second step algorithm\n Global Joint optimization")


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
        help="*.init file to be processed",
        default=None
    )
    return parser.parse_args()



def get_files_with_extension(extension=".init"):
    """
    Retrieves all files in the current directory with the specified extension.
    """
    return sorted([file for file in os.listdir(".") if file.endswith(extension)])

def main():
    """Main function for the script."""

    # Parse arguments
    args = parse_arguments()

    # Read configuration file
    if args.config_file ==None:
        print("Error! give me a cactus config_file")
    else :
        cactus_args = RC.read_config_file(args.config_file)



    # preselect file if specified
    if args.file is not None:

        if not args.file.endswith(".init"):
            print(f"Error!!! the file: {args.file} should be a *.init file")
            print("Exiting...")
            exit()
        files = [args.file]
        print("Processing file: " , files)
    else:
        # Retrieve and process files
        files = RC.get_files_with_extension(".init")
        print(f"Discovered files: {files}")


    for i, file_name in enumerate(files):
        # Skip experiments as specified

        print("#" * 20)
        cactus_args.outfile = file_name

        folder= file_name.split(".")[0]
        profiler = "profiler_ovlp.log"

        current_tolerance = float(RC.read_last_line(os.path.join(folder, profiler)))

        # Check optimization tolerance
        if current_tolerance < cactus_args.optimization_tolerance:
            print("Already optimized.")
            continue

        # Construct and execute the command
        command = (
            f"{cactus_paths.global_optimizer} {cactus_args.outfile} "
            f"{cactus_args.alpha} -1 {cactus_args.max_iterations}"
        )
        os.system(command)
        print("#" * 20)

if __name__ == "__main__":
    main()
