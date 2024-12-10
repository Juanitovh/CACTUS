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
import sparse

###
import numpy as np 
import matplotlib.pyplot as plt
from rich import print
from rich.progress import track
from rich.console import Console
from icecream import ic
###

from icecream import ic 
ic.disable()


from libraries import ascii_art 

from numba import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings

import tqdm
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

from libraries import RotatingStrand as RS
from libraries import WattsonFunction as WF
from libraries import ReadWriteStrands as RWS

from libraries.CactusPaths import cactus_paths

import bz2
import pickle
import _pickle as cPickle

import argparse
parser = argparse.ArgumentParser(description = 'Parameters')
parser.add_argument("-file", type=str, help="list of strands" , default='') 
parser.add_argument("-grid_size", type=float, help="step of metaball space")
parser.add_argument("-iterations", type=int, help="iterations heat" , default = 5)
parser.add_argument("-batch_id", type=int, help="current batch working on" , default=0)
parser.add_argument("-n_batchs", type=int, help="number of batches to split the runn" , default=1)
parser.add_argument("-strand_id", type=str, help="strand id from the separated by commas eg 2,3,11" , default='-1')
parser.add_argument("-missing_axons", type=str, help="list of missing axons to be proccessed" , default=0)
parser.add_argument("-missing_axon_file", type=str, help="list of missing axons to be proccessed" , default="error")
parser.add_argument("-which_bake", type=str, help="file containing the missing strands" , default=None)
parser.add_argument("-q", "--quiet", action="store_true")
parser.add_argument("-f", "--force" , action="store_true",  help="force rewriting")

parser.add_argument("-v", "--verbosity", action="count", default=0)

args = parser.parse_args()



file  = args.file
global_box_step = args.grid_size
global_iterations = args.iterations

part_number = args.batch_id
total_parts =args.n_batchs






@jit(nopython = True)
def sample_spherical(npoints, ndim=3):
    vec = np.random.randn( npoints ,ndim)
    for i in range (npoints):
        vec[i] = vec[i]/np.linalg.norm(vec[i])
    return vec

@jit(nopython= True)
def is_inside_box(p, box):
    v1 = box[0] <= p 
    v2 =  p <= box[1]
    if np.all(v1) and np.all(v2): 
        return True
    return False


@jit(nopython= True)
def get_8_vertices_cube( p0, step , n_neighbour ):
    next_vertex = [0,0,0]
    next_vertex[0] = p0[0] +  ( n_neighbour%2)     *step
    next_vertex[1] = p0[1] +  ( (n_neighbour//2)%2) *step
    next_vertex[2] = p0[2] +  ( (n_neighbour//4)%2) *step
    return np.array(next_vertex)

@jit(nopython= True)
def cube_intersection(p0, step_p0 , box):
    for nn in range(8): 
        next_vertex =  get_8_vertices_cube(p0, step_p0 , nn)
        #if np.all( box[0] <= next_vertex ) and np.all( next_vertex <= box[1]):
        if is_inside_box(next_vertex , box):
            return True
    return False

@jit(nopython= True)
def touches_box(p, box , max_rad ):
    #return cube_intersection( [p[i] - p[3] for i in range(3)] , p[3]*2 , box)
    return cube_intersection( [p[i] - max_rad*2 for i in range(3)] , max_rad*4 , box)


@jit(nopython= True)
def touches_box2(p, box  ):
    #print("current box " , box)
    #return cube_intersection( [p[i] - p[3] for i in range(3)] , p[3]*2 , box)
    return is_inside_box(p[0:3] , box)

### finds the rotation if the vector defining start and end of the strand
def rotation_one_axis(strand):
    e_x = np.array([1,0,0])
    v1 = (strand[0]- strand[ -1])[0:3]
    # print("axis " , v1)
    
    if np.sum(v1) ==0:
        mat = np.eye(3) 
        #print("eye matrix, ball processing" , mat)
    else:
        mat = RS.rotation_matrix_from_vectors(v1, e_x)
    return mat


def get_bb_strand_points(strand ):
    xs ,ys ,zs = strand.T

    extra_rad = 1
    bounding_box = [ np.array([ np.min(xs) , np.min(ys),np.min(zs) ]) - extra_rad*1.5,   np.array([ np.max(xs) , np.max(ys),np.max(zs) ]) + extra_rad*1.5   ]

    return np.array(bounding_box)

def get_vol_box(strand , print_info = False):
    b = get_bb_strand_points(strand )
    if print_info:
        print("bounding box lenght" , b[1]-b[0] , " volume " , np.prod(b[1]-b[0]))

    return np.prod(b[1] - b[0])

def rotation_around_x(theta):
    return np.array([[1, 0, 0],
                     [0, np.cos(theta), -np.sin(theta)],
                     [0, np.sin(theta), np.cos(theta)]])

###
def get_min_rotation_plane(points):
    """
    rotation around x axis because rotate_strand IN LIBRARY works that way with x_axis
    """


    mat_0  = rotation_one_axis(points)
    thetas = np.linspace(0, 2*np.pi, 200)
    volumes = []
    for ti in thetas:
        current_rot = rotation_around_x(ti)
        current_rot = current_rot@mat_0
        pp2 = np.dot(current_rot, points.T).T
        vol = get_vol_box(pp2)
        # print(ti , vol)
        volumes.append(vol)

    mini = np.argmin(volumes)
    good_mat = rotation_around_x(thetas[mini])
    good_mat = good_mat@mat_0

    ppi = np.dot(good_mat, points.T).T

    # print("lenghts box")
    # vol_original = get_vol_box(points , print_info = True)
    vol_original = get_vol_box(points , print_info = False)
    # vol = get_vol_box(ppi, print_info = True)


    #print("#$#$# fraction origin volume ", volumes[mini]/vol_original)
    #print("#$#$# fraction after plane volume ", volumes[mini]/volumes[0] )
    """ print("volumes " , vol_original , volumes[0] ,  volumes[mini] ) """
    return good_mat

def rotate_strands(strands , n_i ):
    #strands = np.copy(strands2)

    mat = get_min_rotation_plane(strands[n_i][: , 0:3])
    for i in range(len(strands)):
        points = strands[i][:, 0:3]
        strands[i][: , 0:3] = mat.dot(points.T).T 
    return mat 

def rotate_strands_oldV(strands , n_i ):
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

@jit(nopython= True , fastmath=True)
def my_norm(x):
    return np.sqrt(np.sum(x**2))

@jit(nopython = True)
def get_sampling_rad(rad , step):
    return np.ceil(4*np.pi * ( (rad/step)**2 )/64)
    #return np.ceil(4*np.pi * ( (rad/step)**2 )/16)
    # return np.ceil(4*np.pi * ( (rad/step)**2 )/4)

@jit(nopython = True)
def fill_grid(grid_temperature ,grid_class ,strands  , bounding_box , new_box , box_step, min_rad , max_rad):
    #start = time.time()
    interpolate_ball = 3
    #extended_radii_box = [bounding_box[0] - max_rad*2 , bounding_box[1] + max_rad*2]
    ex0 = [ v - 2*max_rad for v in bounding_box[0]]
    ex1 = [ v + 2*max_rad for v in bounding_box[1]]
    extended_radii_box = [ex0 , ex1]
    extended_radii_box  = np.array(extended_radii_box).reshape(2,3)
    
    for i in range(len(strands)):
        """ if i%50 == 0: """
        """     print(i , "/", len(strands)) """
        strand = strands[i]
        for j in range(len(strand)):
            ball1 = strand[j]
            if len(strand) == 1:
                ball2 = strand[j]
            elif j+1 == len(strand):
                continue
            else:
                ball2 = strand[j+1]

            
            ball_distance  = my_norm(ball2[0:3]-ball1[0:3])

            interpolate_ball = int(np.maximum(interpolate_ball , ball_distance//(ball1[3] ) ))
            # print("interpolate_ball ", interpolate_ball , "before it was 3")

            step_ball = (ball2-ball1)/interpolate_ball
            for k in range(interpolate_ball+1):
                cp = ball1 + step_ball* k 
                if not touches_box2(cp , extended_radii_box ):
                    continue

                my_hash = np.uint32(hash((i , j ,k)))
                np.random.seed(my_hash)

                #print("I MADE IT")
                current_rad = cp[-1]
                rad_percentage = np.arange(.1*current_rad , 0.95*current_rad , box_step*2)
                #rad_percentage = np.arange(.1*current_rad , .95*current_rad , box_step)
                # rad_percentage = np.arange(.1*current_rad , .95*current_rad , box_step*2)
                if len(rad_percentage) <=2:
                    rad_percentage = np.linspace(.1*current_rad, 0.95*current_rad, 4)

                for per in rad_percentage:
                    surface_sampling = get_sampling_rad(per , box_step)
                    points_sphere = sample_spherical( surface_sampling )*per + (cp[0:3] - bounding_box[0])

                    #if not touches_box(cp , bounding_box):
                    #    continue


                    for p in points_sphere:
                        if is_inside_box(p , new_box ) == False:
                            continue
                        # print(" all is good")
                        p = p/box_step
                        x,y,z = int(p[0]) ,int(p[1]) ,int(p[2])
                        
                        grid_temperature[x][y][z] += 1000
                        grid_class[x][y][z] =  i+1 #1 if i== selected_strand else 2
 


@jit(nopython=True  , nogil = True)
def get_neighbours_3D_index( x,y, z, index):
    return (x+ index%3 -1 , y + index//3 %3-1, z + index//9%3 -1)


@jit(nopython=True  , nogil = True)
def get_neighbours_3D_index_5( x,y, z, index):
    return (x+ index%5 -1 , y + index//5 %5-1, z + index//25%5 -1)

@jit(nopython=True)
def my_argmax(v):
    maxi = -np.inf
    max_index = -1
    for i in range(len(v)):
        if v[i] > maxi:
            maxi = v[i]
            max_index = i
    return max_index 

@jit(nopython=True)
def handle_temperature_old(T , typeT , i,j,k):
    #my_neigh = [get_neighbours_3D_index(i,j,k,index) for index in range(27)]
    current_type = typeT[i,j,k]
    if current_type != 0:
        suma_voxel =0
        #for x,y,z in my_neighbours:
        for index in range(27):
            v = get_neighbours_3D_index(i,j,k, index)
            #x,y,z =my_neigh[index]
            suma_voxel+= T[v] if typeT[v] == current_type else 0
                
        return suma_voxel/27 , current_type

#    types = []
    my_dicto ={}
    for index in range(27):
        v = get_neighbours_3D_index(i,j,k, index)
        inner_type = typeT[v]
        if  inner_type != 0:
            my_dicto[inner_type] =1
        
    types = list(my_dicto.keys())
    #types = list(set(types))
    if len(types)==0 or len(types) >=2:
        return 0,0
    
    dominant_voxel = types[0]
  
    suma_voxel = 0
    for index in range(27):
        v = get_neighbours_3D_index(i,j,k, index)
        if typeT[v] == dominant_voxel:
            suma_voxel+= T[v ] 
    return suma_voxel/27 , dominant_voxel

""" @jit(nopython=True) """
""" def handle_temperature(T , typeT , i,j,k): """
"""     #my_neigh = [get_neighbours_3D_index(i,j,k,index) for index in range(27)] """
"""     current_type = typeT[i,j,k] """
"""     my_dicto ={} """
"""     # my_dicto = defaultdict(lambda:0) """
"""     for index in range(27): """
"""         v = get_neighbours_3D_index(i,j,k, index) """
"""         inner_type = typeT[v] """
"""         if  inner_type != 0: """
"""             try : """
"""                 my_dicto[inner_type] +=1 """
"""             except: """
"""                 my_dicto[inner_type] =1 """
"""          """
"""     types = list(my_dicto.keys()) """
""""""
"""     if len(types)==0 : #or len(types) >=2: """
"""         return 0,0 """
""""""
"""     ii = my_argmax(list(my_dicto.values())) """
"""     predominant = types[ii] """
""""""
"""     dominant_voxel = predominant if predominant > current_type else current_type """
"""    """
"""     suma_voxel = 0 """
"""     for index in range(27): """
"""         v = get_neighbours_3D_index(i,j,k, index) """
"""         if typeT[v] == dominant_voxel: """
"""             suma_voxel+= T[v ]  """
"""     return suma_voxel/27 , dominant_voxel """

@jit(nopython=True)
def my_max(v):
    maxi = -np.inf
    for vi in v:
        if vi > maxi:
            maxi = vi
    return maxi


@jit(nopython=True)
def handle_temperature(T , typeT , i,j,k):
    #my_neigh = [get_neighbours_3D_index(i,j,k,index) for index in range(27)]
    current_type = typeT[i,j,k]
    my_dicto ={}
    # my_dicto = defaultdict(lambda:0)
    for index in range(27):
        v = get_neighbours_3D_index(i,j,k, index)
        inner_type = typeT[v]
        if  inner_type != 0:
            try :
                my_dicto[inner_type] +=1
            except:
                my_dicto[inner_type] =1
        
    types = list(my_dicto.keys())

    if len(types)==0 : #or len(types) >=2:
        return 0,0

    """ if len(types) == 1 and types[0] == current_type: """
    """     return T[i,j,k]*2 , current_type """


    """ ii = my_argmax(list(my_dicto.values())) """
    """ predominant = types[ii] """
    predominant = my_max(types)

    
    if current_type ==0:
        dominant_voxel = predominant
    else:
        dominant_voxel = predominant if  my_dicto[predominant] > 3  else current_type
    """ dominant_voxel = predominant if  my_dicto[predominant] > (my_dicto[current_type] + 2)  else current_type """
    """ dominant_voxel = predominant if (  (predominant > current_type)  and (my_dicto[predominant] > (my_dicto[current_type] + 2)) ) else current_type """

  
    suma_voxel = 0
    for index in range(27):
        v = get_neighbours_3D_index(i,j,k, index)
        if typeT[v] == dominant_voxel:
            suma_voxel+= T[v] 
    return suma_voxel/27 , dominant_voxel


@jit(nopython=True)
def handle_temperature_new(T , typeT , i,j,k):
    #my_neigh = [get_neighbours_3D_index(i,j,k,index) for index in range(27)]
    current_type = typeT[i,j,k]
    my_dicto ={}
    # my_dicto = defaultdict(lambda:0)
    for index in range(27):
        v = get_neighbours_3D_index(i,j,k, index)
        inner_type = typeT[v]
        if  inner_type != 0:
            try :
                my_dicto[inner_type] += T[v]
            except:
                my_dicto[inner_type]  = T[v]
        
    keys    = list(my_dicto.keys())
    my_vals = list(my_dicto.values())

    if len(my_vals)==0 : #or len(types) >=2:
        return 0,0

    """ if len(types) == 1 and types[0] == current_type: """
    """     return T[i,j,k]*2 , current_type """


    predominant = my_argmax(my_vals)
    predominant = keys[predominant]

    """ predominant = types[ii] """
    """ predominant = my_max(types) """

    
    if current_type ==0:
        dominant_voxel = predominant
    else:
        dominant_voxel = predominant if  my_dicto[predominant] > my_dicto[current_type]*1.5 else current_type

  
    """ suma_voxel = 0 """
    """ for index in range(27): """
    """     v = get_neighbours_3D_index(i,j,k, index) """
    """     if typeT[v] == dominant_voxel: """
    """         suma_voxel+= T[v]  """

    return my_dicto[dominant_voxel]/27 , dominant_voxel


from collections import defaultdict
  
      

@jit(nopython=True)
def handle_temperature_cleaning(T , typeT , i,j,k):
    #my_neigh = [get_neighbours_3D_index(i,j,k,index) for index in range(27)]
    if typeT[i,j,k] == 0:
        return  0,0

#    types = []
    my_dicto ={}

    # my_dicto = defaultdict(lambda: 0)
    for index in range(27):
        v = get_neighbours_3D_index(i,j,k, index)
        inner_type = typeT[v] 
        if inner_type != 0:
            my_dicto[inner_type] = 1
        
    types = list(my_dicto.keys())
    #types = list(set(types))
    if len(types)==0 or len(types) >=2:
        return 0,0

    return T[i,j,k] , typeT[i,j,k]


@jit(nopython=True )
def propagate_all_heat_2d(i , grid_temperature , grid_class , lenY, lenZ, my_hanle_temperature , kernel_size =3):
    grid_temperature2d =np.zeros((lenY , lenZ) )
    grid_class2d =np.zeros((lenY , lenZ))
    for j in range(kernel_size//2, lenY-kernel_size//2):
        for k in range(kernel_size//2 , lenZ -kernel_size//2):
            current_temperature , current_type = my_hanle_temperature(grid_temperature , grid_class, i,j,k)
            grid_temperature2d[ j,k] = current_temperature
            grid_class2d[j , k] = current_type


    return grid_temperature2d , grid_class2d


def propagate_heat_iterations(maxIter , grid_temperature , grid_class , my_handle_temperature , kernel_size):
    start = time.time()
    lenX ,lenY , lenZ= grid_temperature.shape
    last_slice_temp = np.zeros((lenY , lenZ))
    last_slice_type = np.zeros((lenY , lenZ))
    for iteration in range(0, maxIter):
       # print(iteration)
        gc.collect()
        last_slice_temp *=0
        last_slice_type *=0
        ####
        
        for i in range(kernel_size//2,lenX-kernel_size//2):
            slice_temp , slice_class = propagate_all_heat_2d(i , grid_temperature , grid_class , lenY , lenZ , my_handle_temperature , kernel_size)

            grid_temperature[i-1 ,: ,:] = last_slice_temp
            grid_class[i-1 ,: ,:] = last_slice_type
            last_slice_temp = slice_temp
            last_slice_type = slice_class


    grid_temperature[0:kernel_size//2] =0
    grid_temperature[-kernel_size//2:] =0

    grid_class[0:kernel_size//2] =0
    grid_class[-kernel_size//2:] =0
    #print("tiempo " , time.time() - start)
    return grid_temperature , grid_class



def get_min_max_rads(strands , selected_strand):
    strands2  = [strands[selected_strand]]
    rads =  [np.max(v[:,3]) for v in strands]
    #rads = strands[selected_strand][: , 3]
    max_rad = np.max(rads)
    min_rad = np.min(rads)
    return min_rad , max_rad


def get_bb_strand(selected_strand , strands , extra_rad , global_bb):
    xs ,ys ,zs = strands[selected_strand][: , 0:3].T

    scale_factor = 2

    bounding_box = [ np.array([ np.min(xs) , np.min(ys),np.min(zs) ]) - extra_rad*scale_factor,   np.array([ np.max(xs) , np.max(ys),np.max(zs) ]) + extra_rad*scale_factor   ]
 
    #for i in range(3):
    #    bounding_box[0][i] = np.maximum(bounding_box[0][i] , global_bb[0][i])
    return np.array(bounding_box)
    #    bounding_box[1][i] = np.minimum(bounding_box[1][i] , global_bb[1][i])




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
        n+=1
    kernel = np.zeros((n,n,n))
    centre = np.array([ n//2 , n//2 , n//2 ])
    rad = n//2+1
    if rad <2:
        rad +=1

    for i in range(n):
        for j in range(n):
            for k in range(n):
                kernel[i,j,k] = in_circle((i,j,k) , centre ,rad)
    return kernel


#%%


import nibabel as nib

def get_mesh_strand(selected_strand  ):
    print("loading mesh strand " , selected_strand)
    zeros_fill=5
    dict_strand = {}

    my_dicto, my_pickle = RWS.get_file_names_pickles_whole(selected_strand)


    # check if the my_pickle file exists
    # print( os.path.isfile(my_pickle) , os.path.isfile(my_dicto) , " dicto and pickle")
    # check if the my_pickle file exists

    if os.path.isfile(my_pickle) and os.path.isfile(my_dicto):
        print("--- "  ,my_pickle , " exists ---")
        return  


    maxIter = global_iterations

    np.random.seed(0)
    color_strand = [np.random.randint(100,256) for i in range(3)]
    dict_strand["color_strand"] = color_strand

    strands , final_lenght_box = RWS.read_generic_list(cactus_paths.optimized_final)
    box_lenght_simulations = final_lenght_box +10

    bb_volume  =  np.array([ (-final_lenght_box/2 , final_lenght_box/2) for i in range(3)]).T
    bb_simulation =  np.array([ (- box_lenght_simulations/2 , box_lenght_simulations/2) for i in range(3)]).T
    dict_strand["bb_volume"] = bb_volume

    rotation_strandi  =  rotate_strands(strands , selected_strand)
    dict_strand["rotation_mat"] = rotation_strandi

    current_rad = np.mean(strands[selected_strand][: , 3])
    dict_strand["current_rad"] = current_rad
## define the working box
    min_rad , max_rad = get_min_max_rads(strands , selected_strand)
    """ max_rad = max_rad*3 """
    
    max_rad = current_rad
    box_step = global_box_step*min_rad
    
        # box_step = 0.1
    dict_strand["box_step"] = box_step


    max_gap = 2.35*maxIter*box_step + max_rad


    ### for each segment of the axon 

    bounding_box = get_bb_strand(selected_strand, strands , max_gap , bb_simulation)
    dict_strand["bounding_box"] = bounding_box

    

    ## for centering in the origin
    new_box = bounding_box - bounding_box[0]

    ##step relative to the minimum axons size
    ic(global_box_step , min_rad)
    #define box dimensions
    box_dim = np.ceil((bounding_box[1] - bounding_box[0] ) / (box_step)).astype(int)


    ## declare the the matrix temperature and classes
    try:
        grid_temperature = np.zeros(box_dim)
        grid_class = np.zeros(box_dim , )
    except : 
        return None

    lenX , lenY , lenZ = grid_class.shape
    

    memory = sys.getsizeof(grid_temperature) + sys.getsizeof(grid_class)
    memory = memory/1e9
    print(f" Processing Memory used  {np.round(memory, 3)} GB   for strand  {selected_strand} ,  || boxstep = {np.round(box_step*1000)} nanoMeter  " )
    # if memory > 1:
    #     print("Memory used " ,  memory , " GB" , " for strand " , selected_strand , "boxstep " , box_step)
    # return
    #print("filling strand: ", selected_strand)
    ic(bounding_box)
    ic(bounding_box[0])

    fill_grid(grid_temperature ,grid_class , strands , bounding_box , new_box , box_step , min_rad , max_rad)

    #save as nii image grid_class

    #nib.Nifti1Image(grid_class, np.eye(4)).to_filename("0.nii.gz")


    

    #save_nifty( selected_strand , "0begin" , grid_class , )
    #print("filled")

    np.random.seed((os.getpid() * int(time.time())) % 123456789)
    color_strand = [np.random.randint(100,256) for i in range(3)]
    dict_strand["color_strand"] = color_strand

    
    dict_strand["maxIter"] = maxIter

    # grid_temperature , grid_class = propagate_heat_iterations(1       , grid_temperature , grid_class , handle_temperature_cleaning , kernel_size=3)
    #save_nifty( selected_strand , "1clean" , grid_class , )
    print("propagating.... " , "|| ", selected_strand)
    grid_temperature , grid_class = propagate_heat_iterations(maxIter , grid_temperature , grid_class , handle_temperature , kernel_size=3)
    print("propagated " , "|| ", selected_strand)

    #nib.Nifti1Image(grid_class, np.eye(4)).to_filename("1.nii.gz")
    # grid_temperature , grid_class = propagate_heat_iterations(1       , grid_temperature , grid_class , handle_temperature_cleaning , kernel_size=3)

    pivot = 1 #np.mean(flat) *.15
    # pivot = .1 #np.mean(flat) *.15

    grid_class[grid_class != (selected_strand+1)] = 0
    
    #empty_places= np.where(grid_class == 0)
    #dict_strand["empty_places"] = empty_places
    #nib.Nifti1Image(grid_class, np.eye(4)).to_filename("11.nii.gz")

    mask_sims = grid_temperature > pivot

    #nib.Nifti1Image(grid_class, np.eye(4)).to_filename("2.nii.gz")
    

    del grid_temperature

    mask_sims = grid_class*mask_sims
    del grid_class

    kernel_7 = circular_kernel(5) # np.ones((7,7,7),np.uint8)
    mask_sims = my_opening(mask_sims , kernel_7 , iterations = 1 )
    mask_sims = my_closing(mask_sims , kernel_7 , iterations = 1 ).astype(bool)

    memory1 = sys.getsizeof(mask_sims)/1e9
    mask_sims = sparse.COO(mask_sims)
    memory2 = sys.getsizeof(mask_sims)/1e9

    print(f"Memory bool   {np.round(memory1,3)} GB ,  Memory compreses: " ,  np.round(memory2,3))

    #mask_vol  =clear_outliers_matrix(mask_sims, bounding_box , bb_volume , box_step , rotation_strandi)
    

    # save pickle dict_strand
    
    print("saving pickle strand " , selected_strand)
    RWS.compressed_pickle_whole(my_dicto, my_pickle , [dict_strand ,mask_sims])
    print("----- done pickle strand " , selected_strand)





#%%

def binary_paste(x, local_itera):
    ## x must be a list of n elements containint two-ples of (local itera dim)
    n = len(x)
    x2 = []
    for i in progressbar.progressbar(range(n//2)):
        x2.append( [x[2*i][kk] + x[2*i+1][kk] for kk in range(local_itera)  ]   )
    if n%2 ==1:
        x2.append(  x[-1])
    return x2



import multiprocessing
from multiprocessing import get_context

prefix_file = file.split('.')[0]

if __name__ =='__main__':
    if args.strand_id == '-1':
        strands_id = [-1]
    else:
        strands_id = list(map(int , args.strand_id.split(',')))
        np.savetxt("last_strands_baked.txt" , strands_id , fmt='%d')


    print(os.getcwd())
    if args.missing_axons ==1:
        strands_id = list(np.loadtxt("0_missing_strands.txt" , dtype=int))

    if args.missing_axon_file == "error":
        ascii_art.print_error_message("Asked to run missing, but there are not any more missing... \n or maybe a problem with the missing file")
        exit()
    elif args.missing_axon_file == "all": 
        ascii_art.print_message("Processing ALL strands togheter")
        strands_id = [-1]
    else:
        ascii_art.print_message("Processing missing strands")
        strands_id = list(np.loadtxt(args.missing_axon_file , dtype=int))


    ascii_art.print_message("Starting the mesh generation")
    #print directory
    print(os.getcwd())

    folder = file.split('.')[0]

    if not os.path.isfile(cactus_paths.optimized_final):
        ascii_art.print_error_message:(f"we need the last version of the optimized file {cactus_paths.optimized_final} "  )


    try:
        os.mkdir(f"meshes/")
        print(f"Creating meshes/ directory")
    except:
        a=1
    try:
        os.mkdir(f"meshes/pickles")
        print(f"Creating meshes/pickle directory")
    except:
        a=1



    ss , final_lenght_box = RWS.read_generic_list(cactus_paths.optimized_final)
    n=len(ss)

    box_lenght_simulations = final_lenght_box +10

    bb_volume  =  np.array([ (-final_lenght_box/2 , final_lenght_box/2) for i in range(3)]).T
    bb_simulation =  np.array([ (- box_lenght_simulations/2 , box_lenght_simulations/2) for i in range(3)]).T


    n_pool = multiprocessing.cpu_count()/2
    n_pool = int(n_pool)


    part_indexes = np.linspace(0,n, total_parts+1).astype(int)
    period = 70

    inputs =  [ i for i in range(n) ]

    if total_parts !=1:
        print(f"subsampling: {args.batch_id}/{args.n_batchs}")
        inputs = inputs[part_indexes[part_number]:part_indexes[part_number+1]]
        inputs= list(reversed(inputs))


    if strands_id[0] != -1:
        inputs = strands_id


        
    if len(inputs ) < n_pool: 
        n_pool = len(inputs)
    print(len(inputs))


    # for ii in inputs:
    #     get_mesh_strand(ii)

    # with get_context("spawn").Pool(n_pool) as pool_obj :

    # pool_obj =  get_context("spawn").Pool(n_pool)
    # res = pool_obj.map(get_mesh_strand, inputs)

    pool_obj =  get_context("spawn").Pool(n_pool) 
    for result in tqdm.tqdm(pool_obj.map(get_mesh_strand, inputs), total=len(inputs)):
        aa=1

    pool_obj.close()
    pool_obj.join()



