#!/usr/bin/env python

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from itertools import count
import copy

data_old = [np.array([1])]
def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Parameters")
    parser.add_argument(
        "-folder", 
        type=str, 
        help="Folder containing log files", 
        default=''
    )
    return parser.parse_args()

def get_best_shape(n):
    """
    Calculate the best layout for subplots based on the number of elements.
    Returns a tuple (rows, columns).
    """
    a = int(np.sqrt(n))
    b = a
    while a * b < n:
        b += 1
    print(f"Best subplot shape: {a} rows, {b} columns")
    return a, b

def find_log_files(folder):
    """
    Find and return a sorted list of .log files in the specified folder.
    """
    return sorted([os.path.join(folder, f) for f in os.listdir(folder) if f.endswith('.log')])

def read_all_logs(files):
    """
    Read data from a list of log files and return it as a list of arrays.
    """
    all_data = []
    for file in files:
        data = np.loadtxt(file)
        if len(data.shape) > 1:
            for i in range(data.shape[1]):
                all_data.append(data[:, i])
        else:
            all_data.append(data)
    return all_data

def animate(i, files, best_shapes, titles):
    """
    Animation function called at each interval.
    Updates the plots if new data is available.
    """
    global data_old
    data = read_all_logs(files)
    

    if data[-1].shape != data_old[-1].shape:
        data_old = copy.copy(data)  # Update the data cache
        plt.cla()  # Clear the current plot
        plt.clf()
        
        print(f"Reading iteration {i}")
        for idx, data_channel in enumerate(data):
            plt.subplot(best_shapes[0], best_shapes[1], idx + 1)
            plt.plot(data_channel, '-^', label='Curent cost')
            
            if idx == len(data) - 1:
                plt.yscale("log")  # Set logarithmic scale for the last subplot
            
            plt.xlabel('Iteration')
            plt.ylabel('Cost')
            plt.title(titles[idx] if idx < len(titles) else f"Plot {idx + 1}")
            plt.legend(loc='upper left')
            plt.tight_layout()
    else:
        print("Data is the same as previous iteration.")
        print("shape of data is ", data[-1].shape)
        print("shape of data_old is ", data_old[-1].shape)

def main():
    """Main function."""
    args = parse_arguments()
    folder = args.folder
    
    # Find all .log files in the folder
    files = find_log_files(folder)
    if not files:
        print("No log files found in the specified folder.")
        return

    print(f"Found files: {files}")
    
    # Determine the best subplot shape
    #best_shapes = get_best_shape(len(files))
    best_shapes = (1,3)
    
    # Read initial data
    data = read_all_logs(files)
    print(f"Loaded data from {len(data)} files.")
    
    
    # Plotting settings
    plt.style.use('fivethirtyeight')
    plt.rcParams["figure.figsize"] = (16, 8)
    
    # Titles for the plots
    titles = ["Length Cost", "Curvature Cost", "Overlapping Cost"]
    data_old = [data[-2]]
    data_old = copy.copy(data_old)
    
    # Start animation
    ani = FuncAnimation(
        plt.gcf(),
        animate,
        fargs=(files, best_shapes, titles),
        interval=2000,
        cache_frame_data=False
    )
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
