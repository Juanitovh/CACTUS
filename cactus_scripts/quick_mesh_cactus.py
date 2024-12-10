#!/usr/bin/env python
#!/usr/bin/python

from libraries import RotatingStrand as RS
from libraries import WattsonFunction as WF
from libraries import ReadWriteStrands as RWS

def tetra_volume(p1, p2 ,p3):
    v321 = p3[0]*p2[1]*p1[2]
    v231 = p2[0]*p3[1]*p1[2]
    v312 = p3[0]*p1[1]*p2[2]
    v132 = p1[0]*p3[1]*p2[2]
    v213 = p2[0]*p1[1]*p3[2]
    v123 = p1[0]*p2[1]*p3[2]
    return (1/6.0)*(-v321 + v231 + v312 - v132 - v213 + v123)

#%%


def get_circle(x,y,z,r,n):
    points = np.zeros((n,3))
    theta = np.linspace(0 , 2*np.pi, n , endpoint= False)
    for i,angle in enumerate(theta):
        points[i,0] = np.cos(angle)*r +x
        points[i,1] = np.sin(angle)*r +y
        points[i,2] = z
    return points.T




def create_circles_from_strand(strand  , n_mesh):
    
    vertices = []
    v1 = strand[0:3,1] - strand[0:3,0]
    vz = np.array([0,0,1])
    points0 = points = get_circle(0,0,0,strand[3,0] , n_mesh)
    mat1 = RS.rotation_matrix_from_vectors(vz, v1)

    circle1 = mat1.dot(points0).T + strand[0:3,0]
    vertices.append(circle1)
    
    for i in range(1,len(strand[0])-1):
        v1 = strand[0:3,i] - strand[0:3,i-1]
        v2 = strand[0:3,i+1] - strand[0:3,i]

        v12 = (v1/np.linalg.norm(v1) + v2/np.linalg.norm(v2))
        
        points0 = get_circle(0,0,0,strand[3,i] , n_mesh)
        mat1 = RS.rotation_matrix_from_vectors(vz, v12)

        
        circlei = mat1.dot(points0).T + strand[0:3 , i]
        vertices.append(circlei)
    
    v1 = strand[0:3,-1] - strand[0:3,-2]
        
    points0  = get_circle(0,0,0,strand[3,-1] , n_mesh)
    
    mat1 = RS.rotation_matrix_from_vectors(vz, v1)

    circle1 = mat1.dot(points0).T + strand[0:3,-1]
    vertices.append(circle1)
    
    
    return vertices
    


def volume_cyl(p1, p2):
    r = (p1[-1] + p2[-1])/2.0
    dista = np.linalg.norm(p1[0:3] - p2[0:3])
    return np.pi*r*r*dista

#%%
import matplotlib.pyplot as plt
import os
import numpy as np
import sys
import pyvista
#os.chdir("/home/juanluis/Documents/strand_optimization/test_cases/jony_list/out/")
#os.chdir("/home/jppl/Desktop/strand_optimization/test_cases/jony_list/out/")
file = sys.argv[1] #"optimized_02950"

try: 
    file_out = sys.argv[2]
except:
    file_out = "test2.ply"

print(file_out)
g_ratio =1
try:
    g_ratio = float(sys.argv[3])
except:
    a=1


""" try: """
"""     #strand_id = int(sys.argv[3]) #2950 """
"""     strand_id = sys.argv[3].split(",") """
"""     strand_id = [int(i) for i in strand_id] """
"""     print(strand_id) """
""" except: """
"""     strand_id = None """

strand_id = None

#strands = np.loadtxt(file)

#%%
##############

strands , box_size = RWS.read_generic_list(file)

max_rads_all = [np.max(v[:,3]) for v in strands]
max_rad = np.max(max_rads_all)
xs = np.concatenate( [v[: , 0] for v in strands])
ys = np.concatenate( [v[: , 1] for v in strands])
zs = np.concatenate( [v[: , 2] for v in strands])

bounding_box = [ np.array([ np.min(xs) , np.min(ys),np.min(zs) ]) - max_rad ,   np.array([ np.max(xs) , np.max(ys),np.max(zs) ]) + max_rad   ]
bounding_box = np.array(bounding_box)
print("--------------------------")
print(bounding_box )

def angle_between(v1, v2):
    #angle between two vectors 3d
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


good_box = []



#%%
##
n_circles = 12
n_circles = 8

n_strands =len(strands)
#n_strands =strands.shape[0]/200


import openmesh as om
import numpy as np
import progressbar

mesh_sphere = None

mesh = om.TriMesh()

rads = []
volumen_approx = 0
if strand_id is None:
    my_range = range(int(n_strands))
else: 
    my_range = strand_id#range(strand_id,strand_id+1)

n_take = 1
""" if len(my_range) > 2000: """
"""     n_take = 3 """
""" elif len(my_range) > 5_000: """
"""     n_take = 5 """
""" elif len(my_range) > 100_000: """
"""     n_take = 100 """
""""""

for i in progressbar.progressbar( my_range):
    #print(i)
    strand = strands[i][::n_take]
    strand[:,3] *= g_ratio
    #strand = strands[i]


    if len(strand) ==1 :
        center = strand[0][0:3]
        rad = strand[0][-1]
        sphere2 = pyvista.Sphere(radius = rad , center = center  , theta_resolution=7 , phi_resolution=7)
        try :
            mesh_sphere += sphere2
        except :
            mesh_sphere = sphere2
        continue

   # strand += np.random.rand(*strand.shape)/20
 #   strand = strands[i*n_control: (i+1)*n_control , :]
  #  rads.append(strand[:,  3]*10 -.1)

    vertices = create_circles_from_strand(strand.T, n_circles)
    #print(vertices)
    
    vers = np.concatenate(vertices)

    list_vert = [ mesh.add_vertex(v) for v in vers]
    n_control_points = len(strand)
    for i in range(n_control_points-1):
        volumen_approx += volume_cyl(strand[i] , strand[i+1])
        for j in range(n_circles):
            v1  = list_vert[n_circles*i + (j)]
            v2  = list_vert[n_circles*i + (j+1)%n_circles]
        
    
            v11  = list_vert[n_circles*(i+1) + (j)]
            v22  = list_vert[n_circles*(i+1) + (j+1)%n_circles]
        
       
            face1 = mesh.add_face(v1,v2,v11)
            face2 = mesh.add_face(v22,v11,v2)
    #        mesh.set_vertex_property('color', v1, [1,0.5,0])
     #       mesh.set_vertex_property('color', v2, [1,0.5,0])
      #      mesh.set_vertex_property('color', v11, [1,0.5,0])
       #     mesh.set_vertex_property('color', v22, [1,0.5,0])
            

    v1 = list_vert[ 0 ]
    for i in range(1,n_circles-1):
        v2 = list_vert[ i ]
        v3 = list_vert[ i+1 ]
        face1 = mesh.add_face(v1,v3,v2)
        
    v1 = list_vert[ -n_circles ]
    
    for i in range(-n_circles+1,-1):
        v2 = list_vert[ i ]
        v3 = list_vert[ i+1 ]
        #face1 = mesh.add_face(v3,v2,v1)
        face1 = mesh.add_face(v3,v1,v2)

            
        



if strand_id is None:
    if mesh_sphere != None: 
        om.write_mesh('test.ply', mesh)

        read_mesh = pyvista.read('test.ply')

        read_mesh = read_mesh + mesh_sphere
        read_mesh.save("test00.ply"  )
        #remove vertex colors

        #dont include normals
        print("saving mesh with spheres")
        read_mesh.save(file_out,  binary = False , recompute_normals = False , texture = None)
    else :

        om.write_mesh(file_out, mesh)
else:
        filename_strand = file_out
        print(filename_strand)
        om.write_mesh(filename_strand, mesh)

## yes
import subprocess
#subprocess.call("meshlab test2.ply" , shell = True)
