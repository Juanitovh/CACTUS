#!/usr/bin/env python

import argparse
import os
import numpy as np
from itertools import product

# Import custom libraries
from libraries import RotatingStrand as RS
from libraries import WattsonFunction as WF
from libraries import ReadWriteStrands as RWS
from libraries import ReadConfigurations as RC
import grid_initialization as InitFace


from libraries import ascii_art 
ascii_art.display_title()

ascii_art.print_message("First step algorithm\n Substrate initialization")


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Cactus Substrate Initialization")
    parser.add_argument(
        "-config_file", 
        type=str, 
        help="Configuration file to initialize several voxel settings. Overrides all parameters.", 
        default=None
    )
    return parser.parse_args()




def main():
    """Main function for creating substrate initialization files."""
    # Display ASCII art



    # Parse arguments
    args = parse_arguments()

    # Read configuration files with morphological features
    # also read a copy to send each of the differnt jobs

    if args.config_file ==None:
        ("Error! give me a cactus config_file")
    else :
        cactus_args = RC.read_config_file(args.config_file, list_format=True)
        cactus_args_case_i = RC.read_config_file(args.config_file, list_format=False)


    # Quick parameter check
    ascii_art.print_message("Parameters read\n from configuration file:")
    print(cactus_args)

    print("")
    input("Press Enter to continue...")

    # Generate all possible combinations of parameters for the different substrates
    parameter_names, iterator_combinations = RC.combination_parameters_iterator(cactus_args)
    number_substrates = len(iterator_combinations)

    ascii_art.print_important_message(f"Proceeding to create \n {number_substrates} substrates")

    # Confirmation for large numbers of substrates
    if number_substrates > 10:
        confirmation = input("Do you want to proceed? (y/n): ").strip().lower()
        if confirmation != "y":
            print("Exiting...")
            exit()

    # Create substrates based on parameter combinations
    for i, combination in enumerate(iterator_combinations):
        # Dynamically update arguments based on the current combination
        for j, key in enumerate(parameter_names):
            setattr(cactus_args_case_i, key, combination[j])

        print("#" * 20)

        # Generate output filename
        cactus_args_case_i.outfile = f"{cactus_args.outfile}_{str(i).zfill(5)}.init"

        # Check if the file already exists
        if os.path.isfile(cactus_args_case_i.outfile):
            print(f"File already exists: {cactus_args_case_i.outfile}")
            continue

        # Start initializing the substrate
        print(f"Creating file: {cactus_args_case_i.outfile}")
        InitFace.create_voxel(cactus_args_case_i)

        print("#" * 20)


if __name__ == "__main__":
    main()
