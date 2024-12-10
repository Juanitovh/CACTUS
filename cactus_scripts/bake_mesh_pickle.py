#!/usr/bin/env python
#!/home/juanluis/anaconda3/bin/python
import gc
import progressbar
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

import bz2
import pickle
import _pickle as cPickle
import sparse
###
import numpy as np 
import matplotlib.pyplot as plt
from rich import print
from rich.progress import track
from rich.console import Console
from icecream import ic
###

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
parser.add_argument("-inn_out", type=str, help="inner outer of both types of strands" , default='outer')
parser.add_argument("-sim_vol", type=str, help="volume or simulation types of mesh" , default='volume')
parser.add_argument("-g_ratio", type=float, help="mean g-ratio used for the axons" , default=0.7)

parser.add_argument("-batch_id", type=int, help="current batch working on" , default=0)
parser.add_argument("-n_batchs", type=int, help="number of batches to split the runn" , default=1)
parser.add_argument("-strand_id", type=str, help="strand id from the separated by commas eg 2,3,11" , default='-1')
parser.add_argument("-which_bake", type=str, help="file containing the missing strands" , default=None)


parser.add_argument("-missing_axons", type=int, help="list of missing axons to be proccessed" , default=0)
parser.add_argument("-missing_axon_file", type=str, help="list of missing axons to be proccessed" , default="error")

parser.add_argument("-n_erode", type=int, help="number of erosions to perfon" , default=0)
parser.add_argument("-colorless", type=int,  help="include color mesh" ,default =1)
parser.add_argument("-q", "--quiet", action="store_true")
parser.add_argument("-f", "--force" , action="store_true",  help="force rewriting")
parser.add_argument("-v", "--verbosity", action="count", default=0)




args = parser.parse_args()

zeros_fill = 5

part_number = args.batch_id
total_parts =args.n_batchs


inn_out = args.inn_out
sim_vol = args.sim_vol

file = args.file

""" try: """
"""     strands_id = list(map(int , args.strand_id.split(','))) """
""" except: """
"""     strands_id =[-1] """


if args.strand_id == '-1':
    strands_id = [-1]
else:
    strands_id = list(map(int , args.strand_id.split(',')))

if args.missing_axon_file is not None:
    if not os.path.isfile(args.missing_axon_file) :
        ascii_art.print_important_message("Asked to run missing, but there are not any more missing... \n or maybe a problem with the missing file")
    else:
        strands_id = list(np.loadtxt(args.missing_axon_file , dtype=int))





max_erosions = int(args.n_erode)



@jit(nopython= True)
def is_inside_box(p, box):
    v1 = box[0] <= p 
    v2 =  p <= box[1]
    if np.all(v1) and np.all(v2): 
        return True
    return False


def rotate_strands(strands , n_i ):
    #strands = np.copy(strands2)

    e_x = np.array([1,0,0])
    v1 = (strands[n_i][0]- strands[n_i][ -1])[0:3]
    
    if np.sum(v1) ==0:
        mat = np.eye(3) 
        #print("eye matrix, ball processing" , mat)
    else:
        mat = RS.rotation_matrix_from_vectors(v1, e_x)

    for i in range(len(strands)):
        points = strands[i][:, 0:3]
        strands[i][: , 0:3] = mat.dot(points.T).T 
        
    return mat 



@jit(nopython = True)
def clear_outliers_matrix(grid_temperature2 , bounding_box , global_bb , step_box , rotation_matrix) :
    grid_temperature = np.copy(grid_temperature2 )

    inv = np.linalg.inv(rotation_matrix)

    a,b,c = grid_temperature.shape
    point = np.zeros(3)
    for i in range(a):
        point[0] = i
        for j in range(b):
            point[1] = j
            for k in range(c):
                point[2] = k
                p1 = point*step_box + bounding_box[0]
                p1 = inv.dot(p1)
                if not is_inside_box(p1, global_bb):
                    grid_temperature[i,j,k]= False

    return grid_temperature 


import cv2
from scipy import ndimage
my_erosion = ndimage.binary_erosion
my_closing = ndimage.binary_closing
my_opening = ndimage.binary_opening

def circular_kernel(n):
    def in_circle(p , centre , rad):
        p = np.array(p)
        dista = np.linalg.norm(p-centre)
        if dista <rad : 
            return True
        return False
    if n%2 == 0:
        n -=1
    kernel = np.zeros((n,n,n))
    centre = np.array([ n//2 , n//2 , n//2 ])
    rad = n//2
    if rad <2:
        rad +=1

    for i in range(n):
        for j in range(n):
            for k in range(n):
                kernel[i,j,k] = in_circle((i,j,k) , centre ,rad)
    return kernel



#%%




        





def get_mesh_strand(selected_strand  ):
        
    ##check if pickle files exist

    my_dicto , my_pickle = RWS.get_file_names_pickles_whole(selected_strand)
    if not os.path.isfile(my_pickle) or not os.path.isfile(my_dicto) :
        print(my_pickle, my_dicto)
        ascii_art.print_error_message(f"pickle info file not found {selected_strand}")
        #missing_strands.append(selected_strand)
        return  selected_strand
    

    file_mesh = RWS.get_file_names_meshes(selected_strand ,  max_erosions , inn_out , sim_vol )

    if os.path.isfile(file_mesh) :
        print("Mesh file found " ,selected_strand)
        return 


    dict_strand, mask3   = RWS.decompress_pickle_whole(my_dicto , my_pickle)
    # mask3 = dict_strand["mask_sims"]  
    mask3 = mask3.todense()

    color_strand = dict_strand["color_strand"]
    bb_volume = dict_strand["bb_volume"]
    rotation_strandi = dict_strand["rotation_mat"]  
    bounding_box = dict_strand["bounding_box"]  
    box_step = dict_strand["box_step"]  
    maxIter = dict_strand["maxIter"]  
    current_rad = dict_strand["current_rad"]


    if max_erosions > 1/box_step -1 :
        print("erosions IS TOO BIG for this grid size")
        sys.exit()

    current_rad = current_rad + box_step*maxIter/2

    def erode_g_ratio_mask(mask, g_ratio):
        rad_kernel_residual = int(current_rad*(1-g_ratio)/box_step*1)*2 ## it has to by x 2 #  1
        
        n_erode_9 = 0
        smaller_kernel = 9

        # print("---- kernels   " , g_ratio  , rad_kernel_residual , current_rad/box_step)
        if rad_kernel_residual > smaller_kernel:
            n_erode_9 = rad_kernel_residual//smaller_kernel
            rad_kernel_residual = rad_kernel_residual - ( n_erode_9*smaller_kernel)

        if rad_kernel_residual <3 and n_erode_9 >=1 :
            n_erode_9 = n_erode_9 -1
            rad_kernel_residual = rad_kernel_residual + smaller_kernel

        print(f"---- g_ratio {np.round(g_ratio,2)}   kernel_1_size  {rad_kernel_residual} -> 1 erosion , kernel_2_size:{smaller_kernel}  ->{n_erode_9} erosions")
        kernel_9 = circular_kernel(smaller_kernel ) ## quick tunning 

        if rad_kernel_residual <3:
            kernel3 = np.ones((3,3,3),np.uint8)
            mask3 = my_erosion(mask , kernel3 , 1)
        else:
            kernel_res = circular_kernel(rad_kernel_residual)
            mask3 = my_erosion(mask , kernel_res , 1)

        if n_erode_9 != 0:
            mask3 = my_erosion(mask3 , kernel_9 , n_erode_9)
            
        return mask3
    #print("n_iterations in " , selected_strand , kernel.shape , n_erode)

    def bake_mesh_and_export(mask3 , name_file):


        vol_file = name_file.split('.')[0] + ".vol"
        f = open(vol_file ,  "w")
        f.write(str(volume_mask) )
        f.close()
        try:
            vertices, faces = measure.marching_cubes(mask3 , step_size = 1)[:2]        

            vertices = vertices*box_step + bounding_box[0]
            my_inv = np.linalg.inv(rotation_strandi)

            vertices = my_inv.dot(vertices.T).T

            faces = np.hstack([ np.concatenate([[3] , v[::-1]]) for v in faces] )
            surf = pv.PolyData(vertices, faces)
            surf = surf.smooth(n_iter=500)

            target_decimate =1- 3/100 #*.5 
            surf = surf.decimate(target_decimate , volume_preservation = True)
            surf = surf.smooth(n_iter=100)


        except :
            print(f"Return empty MESHING {name_file}")
            vertices = np.array([])
            faces = np.array([])
            surf = pv.PolyData([[0,0,0]])

    
        aux2 = RWS.get_file_names_meshes_aux(selected_strand)
        # surf.save(name_file)
        # subprocess.run(f"~/blender/blender --background --python  ~/blender/fix_normals.py -- {name_file} {aux2}"  , shell=True ,capture_output=True)
        # surf = pv.read(aux2)
        # subprocess.run(f"rm {aux2}"  , shell=True ,capture_output=True)

        texture = np.zeros((surf.n_points, 3), np.uint8)

        for coli in range(3):
            texture[:,coli ] = color_strand[coli] 

        if args.colorless == True:
            print("saving colorless")
            surf.save(name_file,  binary = False , recompute_normals = False)
        else: 
            try:
                print("saving rgb")
                surf.point_data['RGB'] = texture
                surf.save(name_file, texture = "RGB" , binary = False , recompute_normals = False)
            except:
                print("error saving color")
                surf.save(name_file,  binary = False , recompute_normals = False)

        return surf

    
    C0 = 0.35 
    C1 = 0.006
    C2 = 0.024

    xx = current_rad
    #log_gratio = (C0+C1*(xx)+C2*np.log(xx)) ## prior taken from novikov paper
    #print("log_gratio " , log_gratio)
    #log_gratio = (C0+C1*(xx)+C2*np.log(xx)) ## prior taken from novikov paper

    #a = 0.24705
    #b = 0.5382

    #my_mean = 1/current_rad *a + b
    my_mean = args.g_ratio

    itera = 0
    np.random.seed(selected_strand)

    current_gratio =np.abs( np.random.normal(my_mean , .04 , 1)[0])
    while (current_gratio >.95 or current_gratio < 0.45) and itera <100:
        current_gratio =np.abs( np.random.normal(my_mean , .04 , 1)[0])
        itera = itera + 1
    if itera > 99:
        current_gratio = my_mean




    if inn_out == "inner":
        mask3 = erode_g_ratio_mask(mask3 , current_gratio)

    kernel_3 = np.ones((3,3,3),np.uint8)
    for i  in range(max_erosions):
        mask3 = my_erosion(mask3 , kernel_3 , 1)


    if sim_vol =='volume':
        mask3  =clear_outliers_matrix(mask3, bounding_box , bb_volume , box_step , rotation_strandi)
        volume_mask = np.count_nonzero(mask3)* (box_step**3)
    else: 
        mask4  =clear_outliers_matrix(mask3, bounding_box , bb_volume , box_step , rotation_strandi)
        volume_mask = np.count_nonzero(mask4)* (box_step**3)
        del mask4

    surf_in  =  bake_mesh_and_export(mask3 , file_mesh)


    return None

#%%






from multiprocessing import get_context
import multiprocessing


file_optimized=  "optimized_final.txt"
if __name__ =='__main__':

    folder = file.split('.')[0]
    print(os.getcwd())

    print(args)

    if args.missing_axon_file == "error":
        ascii_art.print_error_message("Asked to run missing, but there are not any more missing... \n or maybe a problem with the missing file")
        exit()
    elif args.missing_axon_file == "all": 
        ascii_art.print_message("Processing ALL strands togheter")
        strands_id = [-1]
    else:
        ascii_art.print_message("Processing missing strands")
        strands_id = list(np.loadtxt(args.missing_axon_file , dtype=int))


    try:
        os.mkdir(f"meshes")
        print(f"Creating meshes/ directory")
    except:
        print(f"Directory meshes/ already exists")
        a=1

    try:
        os.mkdir(f"meshes/volume")
        print(f"Creating meshes/volume directory")

    except:
        print(f"Directory meshes/volume already exists")
        a=1
    try:
        os.mkdir(f"meshes/simulations")
        print(f"Creating meshes/simulations directory")
    except:
        print(f"Directory meshes/simulations already exists")
        a=1

    try:
        os.mkdir(f"meshes/aux_folder")
        print(f"Creating mesh/aux_folder directory")
    except:
        a=1


    n= RC.count_number_streamlines(cactus_paths.optimized_final)
    print(f"Number of strands {n}")

    n_pool = multiprocessing.cpu_count()//2

    part_indexes = np.linspace(0,n, total_parts+1).astype(int)

    inputs =  [ i for i in range(n) ]


    if total_parts !=1:
        print(f"subsampling: {args.batch_id}/{args.n_batchs}")
        inputs = inputs[part_indexes[part_number]:part_indexes[part_number+1]]


    if strands_id[0] != -1:
        inputs = strands_id


    if len(inputs ) < n_pool: 
        n_pool = len(inputs)
        
    print( "meshing this axons: " , np.array(inputs) , "total " , len(inputs))


    missing_strands = []
    pool_obj = get_context("spawn").Pool(n_pool) 
    for result in tqdm.tqdm(pool_obj.imap_unordered(get_mesh_strand, inputs), total=len(inputs)):
        if result != None:
            missing_strands.append(result)

    print("closing pool")
    pool_obj.close()
    pool_obj.join()


    if len(missing_strands)>0:
        ascii_art.print_error_message(f"Tried to mesh {len(missing_strands)} pickles that are MISSING \n run this same step with flags \'-case missing -substep growth\' ")
    else:
        ascii_art.print_important_message("Lucky you \n everything is fine")





#%%


