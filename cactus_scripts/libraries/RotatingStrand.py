
#import torch
import numpy as np
import math
from numba import jit
import sys


def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """


    # print( np.linalg.norm(vec1) , np.linalg.norm(vec2))
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    if np.isclose(a , b).all() or np.isclose(a,-b).all():
        return np.eye(3)
    # print("a - b " , a , b)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    try:
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    except:
        print(vec1 , vec2)
        sys.exit()
    return rotation_matrix

@jit(nopython= True)
def rotation_matrix_from_vectors_jit(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)

    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

def random_three_vector():
    """
    Generates a random 3D unit vector (direction) with a uniform spherical distribution
    Algo from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    :return:
    """
    phi = np.random.uniform(0,np.pi*2)
    costheta = np.random.uniform(-1,1)

    theta = np.arccos( costheta )
    x = np.sin( theta) * np.cos( phi )
    y = np.sin( theta) * np.sin( phi )
    z = np.cos( theta )
    return (x,y,z)


def curve_torus_az(ax,ay,L, z , N):
    
    offx = np.random.rand()
    offy = np.random.rand()
    
    arr_z = np.linspace(0,2 * np.pi , N)
    
    arrx = ax*(np.cos(arr_z*L + offx) -1)
    arry = ay*np.sin(arr_z*L + offy)
    arr_z = arr_z*z/(2*np.pi)
    return np.array([arrx,arry,arr_z])


def linspace_curved(start , end , ax , ay , L, N):
    
    norm = np.linalg.norm(end- start)
    
    points = curve_torus_az(ax,ay, L , norm , N)
    
    vec1 = np.array([0,0,1])
    vec2 = end - start
    mat = rotation_matrix_from_vectors(vec1, vec2)
    
    points = mat.dot(points)
    
    #points = torch.tensor(points)
    
    points[0 , :]  =  points[0 , :] + start[0]
    points[1 , :]  =  points[1 , :] + start[1]
    points[2 , :]  =  points[2 , :] + start[2]
    
    return points


def linspace_curved_with_rad(start , end , ax , ay , L, N , radius):
    
    norm = np.linalg.norm(end- start)
    
    points = curve_torus_az(ax,ay, L , norm , N)
    
    vec1 = np.array([0,0,1])
    vec2 = end - start
    mat = rotation_matrix_from_vectors(vec1, vec2)
    
    points = mat.dot(points)
    
    #points = torch.tensor(points)
    
    points[0 , :]  =  points[0 , :] + start[0]
    points[1 , :]  =  points[1 , :] + start[1]
    points[2 , :]  =  points[2 , :] + start[2]
    
    points2 = np.zeros((4,N))
    points2[0:3 , : ] = points
    points2[3 , : ] = radius
    
    return points2

def linspace_curved_cover(start , end , N):
    
    return linspace_curved(start , end , 0.1 , 0.1  , 1, N).T



def sample_uniformely_strand(strand , step_size = 10 ):
    """
    give an skeleton point every step_size units to subsample or up-sample.
    equivalent to np.arrange but on a generic n-d curve with last dimension == radius
    """
    points = strand
    n = len(strand)-1
    cum_step = 0
    
    points_skelet = [np.array(points[0])]
    distances = []
    for i in range(n):
        v_dir = points[i+1] - points[i]
        v_dir= v_dir[0:-1]
        norm_v = np.linalg.norm(v_dir)

        if cum_step + norm_v >= step_size:
            v_unit = v_dir/norm_v

            cum_step2 = step_size - cum_step
   #         print("------ " , cum_step , cum_step2 , cum_step2 + cum_step)        
            while  cum_step2 <= norm_v:
                new_pos = points[i][0:-1] + cum_step2 *v_unit
                cum_step2 += step_size 
                new_pp = np.zeros_like(points[i])
                new_pp[0:-1] = new_pos
                new_pp[-1] = points[i][-1]
                points_skelet.append(new_pp)
                distances.append(np.linalg.norm(points_skelet[-1][0:3] -  points_skelet[-2][0:3]))

            cum_step = norm_v - (cum_step2 -step_size)
        else: 
            cum_step += norm_v

    #points_skelet.append(points[-1])
    return np.array(points_skelet) 

