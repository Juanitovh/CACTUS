import numpy as np
###
import numpy as np 
import matplotlib.pyplot as plt
from rich import print
from rich.progress import track
from rich.console import Console
from icecream import ic
###

import bz2
import pickle
import _pickle as cPickle
import sparse
import os 


from libraries import wavy_radii as WR
import sys
## create streamline from start and end points
### creating a streamline with radius from start to end with parameters of the numbres of control points


def create_streamline(pos_start, pos_end , PerExtraDist = 0.4 , wavy_radii = False , lenght_stick = 5 , lenght_ball = 5 , min_rad = .5, max_rad = 4):

    start ,end  ,r = pos_start[0:3] , pos_end[0:3] , pos_start[-1]
    start, end = np.array(start) , np.array(end)

    len_strand = np.linalg.norm(start -end)

    #end += .005
    n_control = np.linalg.norm(start -end)/(1.25*r) * (1 + PerExtraDist)
    n_control = np.linalg.norm(start -end)/(1.25*(r**.5)) * (0.6 + PerExtraDist)
    #n_control = np.linalg.norm(start -end)/(1.25*r) * (1 + 0)
    #n_control = np.linalg.norm(start -end)/(.25*r) * (1 + 0)

    step_distance = np.linalg.norm(start -end)/n_control
    """ ic(step_distance , n_control) """
    if n_control <1:
        return np.array([pos_start]) , 1

    n_control = np.minimum(1300 , n_control)
    n_control = int(n_control)
    start = np.array([*start, r])
    end = np.array([*end, r])
    # cylinder = np.linspace(start , end ,  num = int(n_control))
    cylinder = np.linspace(start , end ,  num = int(n_control) , endpoint = True)

    if wavy_radii ==True: 

        id_capsule, final_radii = WR.oscillating_radii( lenght_strand=len_strand, n_control_points = n_control, lenght_ball=lenght_ball, lenght_stick=lenght_stick, plot=False)
        #sinus_rad = np.arange(len(cylinder))/len(cylinder) *np.linalg.norm(start -end)
        var_val = (min_rad/r)**.05     #logaritmic scale for radii increase in 
        min_val = (min_rad/max_rad)**.05
        max_val = (min_rad/min_rad)**.05


        dynamic_range_rad = var_val  - min_val
        #dynamic_range_rad = dynamic_range_rad / (max_val - min_val)   ## this is number 5 
        """ dynamic_range_rad = dynamic_range_rad / (max_val - min_val) + .2 #.15 # this is number 7 """
        dynamic_range_rad = dynamic_range_rad / (max_val - min_val) + .05 #.15 # this is number 7

        #p = ( 1- ( r - min_rad)/(max_rad - min_rad)  ) *.8
        #np.random.uniform(p-.2 , p )         +
        scale_rad = 1  #  0.1


        #final2 = (final_radii-1) *scale_rad + .9 # now it's in interval (0.75 , 1.5)
        final2 = (final_radii-1) *dynamic_range_rad  # .9 # now it's in interval (0.75 , 1.5)
        final2 = final2 - np.min(final2)*1 +.9
        """ final2 = final2 - np.min(final2)*1.2 + 0. """
        #mean_final = np.mean(final2)*1.5
        #wb =(1-mean_final) 
        """ if wb < 0: """
        """     print('wb is negative') """
        """ else: """
        """     wb= .1 """
        """ final2 = final2 +  wb """
        """ ic(dynamic_range_rad , mean_final , wb ) """



        np.random.seed()
        """ ic(final_radii.shape) """
        """ ic(cylinder.shape) """
        #sinus_rad = np.sin(sinus_rad/(.5*np.pi) + np.random.rand()*np.pi*2 ) * r/3
        """ cylinder[:,3] = final2 * r """
        before_correction = final2 * r

        diff = np.mean(before_correction) - r
        final3 = before_correction - diff

        cylinder[:,3] = final3
        plt.plot(final2*r , 'green' , label = 'radii axon')
        plt.axhline(np.mean(final2*r) , color = 'red' , label = "mean radii axon")
        plt.axhline(r , color = 'black' , label = "original radii")
        plt.plot(final3 , color = 'cyan' , label = 'corrected radii axon')
        plt.axhline(np.mean(final3) , color = 'purple' , label = 'mean corrected ' , linewidth = 5 ,linestyle = '-.')
        plt.legend()
        plt.show()


    #start[0] -= 3
    #end[0] += 3
    # cylinder = np.array([start , *cylinder, end])
    # noise = np.random.normal(0, .1 , size = (cylinder.shape[0]-2 ,3))
    # cylinder[1:-1 , 0:3 ] +=  noise
    return cylinder , n_control



#write generic list python
def write_strand_list_with_box(cylinders , out_file = 'jony_cylist2.txt' , box_size = -1):
    f = open(out_file , 'w')
    f.write(f'{box_size}\n')
    f.write(f'{len(cylinders)}\n')
    for cyli in cylinders:
        f.write(f'{len(cyli)}\n')
        for point in cyli:
            f.write("{} {} {} {}\n".format(*point))
    
    f.close()


def write_strand_list_with_box_jony(cylinders , out_file = 'jony_cylist2.txt' , box_size = -1  , scale_rad = 1):
    f = open(out_file , 'w')
    f.write('0.001\n')
    for cyli in cylinders:
        x,y,z,r = cyli[0]
        x2,y2,z2,r2 = cyli[-1]
        f.write(f'{x} {y} {z} {x2} {y2} {z2} {r2*scale_rad}\n')
    f.close()



def give_float_list(line):
    nums= [float(v) for v in line.split()]
    return nums

def read_generic_list(file):
    f = open(file , 'r')
    box_size = float(f.readline().split()[0])
    n_s = int(f.readline().split()[0])
    lista_cyl = []
    
    for i in range(int(n_s)):
        n_control = int( f.readline().split()[0] )
        cylinder = []
        for j in range(n_control):
            line = f.readline() 
            point = give_float_list(line)
            cylinder.append(point)
        lista_cyl.append(np.array(cylinder))        
    f.close()

    return lista_cyl , box_size


def get_file_names_meshes(selected_strand ,  n_erosion, type_strand , folder):
    zeros_fill = 5
    ff2 = f"meshes/simulations/strand_{str(selected_strand).zfill(zeros_fill)}_{type_strand}_erode_{n_erosion}.ply" 
    return ff2

def get_file_names_meshes_aux(selected_strand ):
    aux2 = f"meshes/aux_folder/aux2_{selected_strand}.ply"
    return aux2



def get_file_names_pickles_whole(selected_strand  ):
    zeros_fill = 5
    my_pickle = f"meshes/pickles/strand_{str(selected_strand).zfill(zeros_fill)}_dict.npz"
    my_dicto = f"meshes/pickles/strand_{str(selected_strand).zfill(zeros_fill)}_dict.pbz2"
    return my_dicto, my_pickle



def compressed_pickle_whole(my_dicto , my_pickle, data):
    #print(my_dicto , my_pickle, data)
    pikd = open(my_dicto , 'wb')
    pickle.dump(data[0], pikd)
    pikd.close()
    sparse.save_npz(my_pickle, data[1])

def decompress_pickle_whole(my_dicto , my_pickle):
    pikd = open(my_dicto, 'rb')
    dicto = pickle.load(pikd)
    pikd.close()
    mat1 = sparse.load_npz(my_pickle)
    return  dicto, mat1


