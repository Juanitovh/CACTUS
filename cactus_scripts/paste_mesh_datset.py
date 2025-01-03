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
###
import numpy as np 
import matplotlib.pyplot as plt
from rich import print
from rich.progress import track
from rich.console import Console
from icecream import ic
###

import subprocess

from numba import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings

import tqdm
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

from libraries import RotatingStrand as RS
from libraries import WattsonFunction as WF
from libraries import ReadWriteStrands as RWS

import argparse
parser = argparse.ArgumentParser(description = 'Parameters')
parser.add_argument("-file", type=str, help="list of strands" , default='') 
parser.add_argument("-n_erode", type=str, help="erosion index")
parser.add_argument("-sim_vol", type=str, help="folder of pasting" , default = "simulations")
parser.add_argument("-inn_out", type=str, help="inner or outer" , default="outer")
parser.add_argument("-decimate", type=float, help="decimation [0,1)" , default=0)
parser.add_argument("-patch_flip", type=int, help="0/1" , default=1)
parser.add_argument("-colorless", action='store_true', help="colorless in vertex")

#parser.add_argument("-q", "--quiet", action="store_true")
#parser.add_argument("-f", "--force" , action="store_true",  help="force rewriting")

#parser.add_argument("-v", "--verbosity", action="count", default=0)


args = parser.parse_args()
print(args)


def getter_filename(selected_strand , inn_out = 'outer' , sim_vol  = 'simulations'  , n_erode = 0):
    zeros_fill =5
    ff2 = f"{sim_vol}/strand_{str(selected_strand).zfill(zeros_fill)}_{inn_out}.ply" 
    return ff2

def get_final_file( inn_out  = 'outer' , sim_vol  = 'simulations'  , n_erode = 0):
    f2 = f"{sim_vol}/hp_{inn_out}.ply" 
    return f2



def paste_half(x):
    ## x must be a list of n elements containint two-ples of (local itera dim)
    n = len(x)
    x2 = []
    for i in progressbar.progressbar(range(n//2)):
        x2.append( x[2*i] + x[2*i+1]  )
    if n%2 ==1:
        x2.append(  x[-1])
    return x2

def paste_binary(x):
    while len(x)>=2:
        x = paste_half(x)
    return x


from numba import jit

@jit(nopython=True)
def tetra_volume2(p1, p2 ,p3):
    v321 = p3[0]*p2[1]*p1[2]
    v231 = p2[0]*p3[1]*p1[2]
    v312 = p3[0]*p1[1]*p2[2]
    v132 = p1[0]*p3[1]*p2[2]
    v213 = p2[0]*p1[1]*p3[2]
    v123 = p1[0]*p2[1]*p3[2]
    return (1/6.0)*(-v321 + v231 + v312 - v132 - v213 + v123)



@jit(nopython=True)
def tetra_volume(p1, p2, p3):
    # Compute the volume of the tetrahedron formed by the origin and the triangle (p1, p2, p3)
    return np.dot(p1, np.cross(p2, p3)) / 6.0


@jit(nopython=True)
def compute_vol(vers , all_faces):

    """ vers = my_surface.points """
    """ all_faces = my_surface.faces.reshape(-1 , 4) """

    volume_tetra = 0
    for i in range(len(all_faces)):
        _,a,b,c = all_faces[i]
        volume_tetra +=  tetra_volume(vers[a] ,  vers[c], vers[b]  )

    return volume_tetra

def compute_true_volume(my_surface):

    vers = my_surface.points
    all_faces = my_surface.faces.reshape(-1 , 4)

    volume_tetra = 0
    for i in progressbar.progressbar(range(len(all_faces))):
        _,a,b,c = all_faces[i]
        volume_tetra +=  tetra_volume(vers[a] ,  vers[c], vers[b]  )


    min_box = bb_volume[0] 
    max_box = bb_volume[1]

    vv = max_box - min_box

    volume = vv[0] * vv[1] *vv[2]
    final_icvf = volume_tetra/volume
    print("ICVF ===== " , final_icvf)
    return final_icvf

#ss , final_lenght_box = RWS.read_generic_list(args.file)

    #read the file first two lines
with open(args.file) as f:
    first_line = f.readline()
    second_line = f.readline()

nn = int(second_line)
final_lenght_box = float(first_line)

ffinal = get_final_file(args.inn_out , args.sim_vol , args.n_erode)
print("------------ Saving final file " , ffinal)

bb_volume  =  np.array([ (- final_lenght_box/2 , final_lenght_box/2) for i in range(3)]).T
import os
print(os.getcwd())
meshes = []
missings = []
bar = progressbar.progressbar
failed_n_verts = []

bad_meshes = []

all_volumes = []

""" nn = 1000 """
for i in track(range(nn)):
    file = getter_filename(i , args.inn_out , args.sim_vol , args.n_erode)
#    print("fixing normals with blender")
#    subprocess.run(f"blender --background --python  /media/juanluis/HDD/strand_optimization/fix_normals.py -- {file} {file}" , shell = True , capture_output=True)
#    print("done")

    if not os.path.isfile(file):
        
        print("error missing strands " , i )
        print(file)
        missings.append(i)
    else:

        try:
            surf = pv.read(file)
        except:
            surf = None
            print("bad mesh" , file)
            bad_meshes.append(file)

# sum up the volumes to get the total volume

        total_volume = compute_vol(np.array(surf.points) , np.array(surf.faces).reshape(-1 , 4) )

        if args.patch_flip and total_volume > 0:
            print("flipping normals" , file)
            surf.flip_normals()
            surf.save(file ,  binary = False , recompute_normals = False )

            total_volume = compute_vol(np.array(surf.points) , np.array(surf.faces).reshape(-1 , 4) )

        all_volumes.append(total_volume)
        #surf4 = surf3.clean()
        #triangulate
        #surf.is_all_triangles = true
        #ic(surf.is_all_triangles)

        #surf  = trimesh.load(aux2)

        if args.decimate != 0:
            try:
                ## args.decimate is target reduction
                surf2 = surf.decimate(args.decimate , volume_preservation = True)
                #surf2 = surf.decimate_pro(args.decimate , )

            except :
                """ print("decimation failed " , file) """
                surf2 = surf
                failed_n_verts.append(len(surf.points))
                pass


            surf = surf2


            color_strand = [np.random.randint(100,256) for i in range(3)]
            texture = np.zeros((surf.n_points, 3), np.uint8)
            for coli in range(3):
                texture[:,coli ] = color_strand[coli] 
            try:
                surf.point_data['RGBA'] = texture
            except:
                surf.point_data['RGB'] = texture


        if len(surf.points) != 0 :
            meshes.append(surf)

ic(len(failed_n_verts))


print("#######################")
print("bad meshes" , len(bad_meshes) , bad_meshes[:10])
print("#######################")


#plt.hist(failed_n_verts)
#plt.show()




print("Missing strands" , len(missings) , missings[:10])

all_volumes = np.array(all_volumes)

min_box = bb_volume[0] 
max_box = bb_volume[1]

vv = max_box - min_box

volume_tetra = np.sum(-all_volumes)
volume = vv[0] * vv[1] *vv[2]
final_icvf = volume_tetra/volume
print("ICVF ===== " , final_icvf)


plt.subplot(121)
#title for both plots
plt.suptitle(f"ESTIMATED ICVF : {np.round(final_icvf*100, 2)} %")
plt.title("Volumes of axons ")
plt.plot(-all_volumes , '.')
plt.ylabel("Volume $\mu m^3$")
plt.xlabel("index fibre")
plt.subplot(122)
plt.title("approx rad")
plt.plot(((-all_volumes)/(final_lenght_box *np.pi) )**.5 , '.')
plt.ylabel("Approx radius $\mu m$")
plt.xlabel("index fibre")

plt.show()

#%%
positive_vols = np.where(np.array(all_volumes) > 0)[0]

ic(positive_vols)

if len(positive_vols) !=0:
    np.savetxt("positive_vols.txt" , positive_vols , fmt = "%d")
    sys.exit()


sys.exit()
#%%

ffinal = get_final_file(args.inn_out , args.sim_vol , args.n_erode)

ic(ffinal)

meshes= paste_binary(meshes)

print("saving mesh")
print(ffinal)

if args.colorless :
    meshes[0].save(ffinal ,  binary = False , recompute_normals = False )
else: 
    try:
        meshes[0].save(ffinal , texture = 'RGBA', binary = False , recompute_normals = False)
    except:
        meshes[0].save(ffinal , texture = 'RGB', binary = False , recompute_normals = False)
print("saving mesh")

print("Pasted")


#for i in range(max_erosions):
#print("Volume erosion " , i )
#compute_true_volume(data_vols_out[0][i])

# volume = compute_true_volume(meshes[0])

# def get_final_file_vol( inn_out = 'outer' , sim_vol  = 'simulations'  , n_erode = 0):
#     f2 = f"meshes_final/{sim_vol}/hp_{sim_vol}_{inn_out}_erode_{n_erode}_triangles.icvf" 
#     return f2
#
# f_vol = get_final_file_vol(args.inn_out , args.sim_vol , args.n_erode)
# np.savetxt(f_vol , np.array([volume]))
#
