#!/usr/bin/env python
#!/home/juanluis/anaconda3/bin/python
import os
import numpy as np


import numpy as np

import argparse


from libraries import ReadWriteStrands as RWS




def get_bb_strand(strand ):
    xs ,ys ,zs , rs = strand.T
    extra_rad = np.max(rs)
    bounding_box = [ np.array([ np.min(xs) , np.min(ys),np.min(zs) ]) - extra_rad*1.5,   np.array([ np.max(xs) , np.max(ys),np.max(zs) ]) + extra_rad*1.5   ]
    return np.array(bounding_box)

def get_bb_all_strands(strands ):
    
    min_corners = []
    max_corners= []

    for strand in strands:
        low_corner , max_corner = get_bb_strand(np.array(strand))
        min_corners.append(low_corner)
        max_corners.append(max_corner)

    min_corners = np.array(min_corners)
    max_corners = np.array(max_corners)

    low_corner = np.min(min_corners,axis=0)
    max_corner = np.max(max_corners,axis=0)

    mid_point = (max_corner + low_corner)/2

    return np.array([low_corner,max_corner]) , mid_point


import progressbar
import matplotlib.pyplot as plt
def read_nfg(my_folder):


    files = sorted(os.listdir(my_folder))

    # files = files[1:10]
# example= "strand_00-01-r0.015000.txt"


    strands = []
    bar = progressbar.ProgressBar(max_value=len(files))
    my_rads = []
    for f in bar(files):

        if not f.endswith(".txt"):
            continue
        
        strand = np.loadtxt(my_folder + f)[::1]
        
        radius = float(f.split(".txt")[0].split('-r')[1]) ###get radius from name
        rads = np.repeat(radius, len(strand)).reshape(len(strand),1)
        my_rads.append(radius)

        strand = np.append(strand, rads , axis=1)
        strands.append(strand*1000)

    index = list(reversed(np.argsort(my_rads)))
    plt.hist(my_rads)
    plt.axvline(np.mean(my_rads) , color = "r")
    # plt.show()

    strands = [np.array(strands[i]) for i in index]

    return strands 


def center_around_origin(strands):

    bb , center = get_bb_all_strands(strands)

    print("old bounding box" , bb)

    for i in range(len(strands)):
        strands[i][:,0:3] -= center[0:3]
        # strands[i][:,3] -= center[3]

    bb , center = get_bb_all_strands(strands)
    print("new bounding box" , bb)
    lenght_side = np.max( bb[1] - bb[0] )
    print("length side" , lenght_side)

    return strands , lenght_side


def process_my_strands(folder , outfile):
    strands = read_nfg(folder)
    strands , lenght_side = center_around_origin(strands)
    RWS.write_strand_list_with_box(strands , out_file = outfile , box_size = lenght_side)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = 'Parameters')
    parser.add_argument("-folder", type=str, help="folder with strands nfg style" , default='') 
    parser.add_argument("-outfile", type=str, help="outfile with list of strands" , default='') 

    args = parser.parse_args()

    process_my_strands(args.folder , args.outfile)




