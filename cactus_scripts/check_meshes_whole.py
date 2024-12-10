#!/usr/bin/env python
#!/home/juanluis/anaconda3/bin/python
import gc
import progressbar
import nibabel as nib
import time
from numba import jit
import os
import numpy as np
import subprocess
from skimage import measure
import pyvista as pv
import sys 
import multiprocessing                                                                                  
import matplotlib.pyplot as plt
import sys

import multiprocessing
import tqdm

from numba import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings

import tqdm
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

from libraries import RotatingStrand as RS
from libraries import WattsonFunction as WF
from libraries import ReadWriteStrands as RWS
from libraries import ReadConfigurations as RC
from libraries import ascii_art
from libraries.CactusPaths import cactus_paths

import argparse
parser = argparse.ArgumentParser(description = 'Parameters')
parser.add_argument("-file", type=str, help="list of strands" , default='') 
parser.add_argument("-folder", type=str, help="folder_analyze" , default='') 
parser.add_argument("-outfile", type=str, help="file to save indexs" , default='0_missing_strands.txt') 
parser.add_argument("-n_erode", type=int, help="erosion number" , default=0) 
parser.add_argument("-inn_out", type=str, help="inner/outer" , default='outer')
parser.add_argument("-delete_pickles",  type=int , help="wheter to delete pickles" , default=0)
parser.add_argument("-missing_axon_file", type=str, help="list of missing axons to be proccessed" , default=None)



args = parser.parse_args()
print(args)



outfile = args.outfile


## read file until find line equal to string "end_header"
def read_ply_until_end_header(file_name):
    n_vertex = -1
    n_faces = -1
    
    total_lines = 0

    with open(file_name, 'r') as f:
        for line in f:
            line2 = line.strip()
            words = line2.split(' ')

            if words[0]=='element':
                if words[1]=='vertex':
                    n_vertex = int(words[2])
                if words[1]=='face':
                    n_faces = int(words[2])


            if line.strip() == "end_header":
                break
            else:
                continue
        count = 0

        for line in f:
            count+=1

    #print("total lines", count , "n_vertex", n_vertex , "n_faces", n_faces , "sum total", n_vertex+n_faces)
    if count != n_vertex + n_faces :
        print("error in reading ply file")
        return -1

    return 1


#count number lines in file



prefix_file = args.file.split('.')[0]

## read first two lines of the file prefix_file


n_strands = RC.count_number_streamlines(cactus_paths.optimized_final)

#strands , final_lenght_box = RWS.read_generic_list(args.file)


#n_strands = len(strands)

def try_read_ply(i):

    strand_file  = RWS.get_file_names_meshes(i , args.n_erode , args.inn_out , args.folder)
    #check if strand_file exists
    idi = None
    if not os.path.isfile(strand_file):
        missing_strands.append(i)
        #print(f"{strand_file} is missing")

        if args.delete_pickles == 1:
            pickle_file , dict_file = get_file_names_pickles_whole(i)
            if os.path.isfile(pickle_file):
                os.remove(pickle_file)
                print(f"pickle {pickle_file} deleted")

        return i
    """ try : """
    """     res = read_ply_until_end_header(strand_file) """
    """     if res ==-1: """
    """         print(f"{strand_file} error ") """
    """         #delete file strand_file  """
    """         #os.delete(strand_file) """
    """"""
    """         return i """
    """ except: """
    """     print(f"unable to read {i}") """
    """     return i """

    return None

missing_strands = []

bar = progressbar.progressbar

n_pool = np.minimum(n_strands , 20)
pool_obj= multiprocessing.Pool(10)

inputs = range(0, 9280)
inputs = range(0, n_strands)
#inputs = range(9580, 9580+50)

for result in tqdm.tqdm(pool_obj.imap_unordered(try_read_ply, inputs), total=len(inputs)):
#for i in bar(inputs):
#    result = try_read_ply(i)
    if result != None:
        missing_strands.append(result)


#pool_obj.close()
#pool_obj.join()

print("pool closed")

if len(missing_strands)>0:
    ascii_art.print_important_message(f"There are {len(missing_strands)} meshes MISSING \n that need to be proccesed ")
else:
    ascii_art.print_important_message("Lucky you \n everything is fine")




if len(missing_strands) ==0 :
    print("eveything is fine")
    if os.path.isfile(outfile):
        os.remove(outfile)
else:
    ascii_art.print_message("Saving files to be proccesed")
    print("missing strands" , missing_strands[0 : np.min([10 , len(missing_strands)])])
    np.savetxt(outfile , missing_strands , fmt='%d')


