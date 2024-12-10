#!/usr/bin/env python

import os
import seaborn as sns

from numba import jit
import multiprocessing                                                                                  

from libraries import WattsonFunction as WF
from libraries import ReadWriteStrands as RWS
from libraries import ReadConfigurations as RC
max_limit = 10
from sklearn.neighbors import BallTree
import time
import progressbar
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
from numba import jit

import numpy as np
import math as m

from sklearn.neighbors import NearestNeighbors
from icecream import ic

#os.chdir('/media/juanluis/HDD/strand_optimization/better_prioirs')
import sys
# sys.argv = ['-q']

import argparse

import psutil
max_cores = psutil.cpu_count(logical = False)//2

parser = argparse.ArgumentParser()
parser.add_argument("-icvf", type=float, help="Intra Cellullar volume fraction" , default =0.50)
parser.add_argument("-outfile", type=str, help="filename output" , default = 'dummy_cross.txt')
parser.add_argument("-outfile2", type=str, help="filename output" , default = 'dummy_cross2.txt')
parser.add_argument("-lenght_side", type=float, help="distance between two faces" , default=50)
parser.add_argument("-dispersion", type=float, help="average dispersion degree" , default=5)
parser.add_argument("-n_bundles" , type=int , help= "number of fibers", default = 1)
parser.add_argument("-crossing_angle" , type=float , help= "crossing_angle fibers", default = 45)
parser.add_argument("-overlap" , type=float , help= "overlap between fibers, percentage", default = .5)
parser.add_argument("-n_cores" , type=int , help= "number of cores", default = max_cores)

parser.add_argument("-gamma_theta", type=float, help="theta parameter in gamma distribution" , default=-1)
parser.add_argument("-gamma_kappa", type=float, help="kappa parameter  in gamma distribution" , default=-1)
#parser.add_argument("-wavy_radii", type=int, help="kappa parameter  in gamma distribution" , default=0)
parser.add_argument("-bias", type=float, help="how much to substract to radii" , default=0.0)

parser.add_argument("-min_rad", type=float, help="min axon radii" , default=.15)
parser.add_argument("-max_rad", type=float, help="min axon radii" , default=4)


parser.add_argument("-v", "--verbose", action="count", default=0)

parser.add_argument("-config_file", type=str, help="filename config file overrides all parameters" , default = None)


args = parser.parse_args()








# distance to work with the ball tree, 
@jit(nopython=True , fastmath= True)
def my_distance(v,v2):
    res =  sqrt((v[0]-v2[0])*(v[0]-v2[0]) + (v[1]-v2[1])*(v[1]-v2[1])) - v[2] - v2[2] + max_limit
    return res
# cosine angle two vectors
def cosine_2vectos(vec_front , vec_skelet):
    vec_front = vec_front/np.linalg.norm(vec_front)
    vec_skelet = vec_skelet/np.linalg.norm(vec_skelet)
    return np.dot(vec_front , vec_skelet)

### get absolute angle between two vectors , angle is between 0 and 90!
def get_angle( vector , reference): 
    angle = cosine_2vectos(vector, reference) 
    my_angle = np.arccos(angle)/np.pi*180   # arccos returns number in [0,pi]
    my_angle = np.minimum(my_angle , np.abs( 180-my_angle))
    return my_angle


def plot_circle(o= [0,0] , rad = 1 , col = None):
    """ theta = np.linspace(0,np.pi*2 , 10) """
    theta = np.linspace(0,np.pi*2 , 30)
    x = [rad * np.cos(ti) + o[0] for ti in theta]
    y = [rad * np.sin(ti) + o[1] for ti in theta]
    
    if col == None:
        plt.plot(x,y, '-')
    else :
        plt.plot(x,y, col+'-' , linewidth = 2)


def plot_circle_array(vec):
    for xi,yi,ri in vec:
        """ plot_circle((xi,yi),ri) """
        plot_circle((xi,yi),ri , 'k')


# check if circle is inside vox
@jit(nopython = True)
def is_inside_box(circle , edges):
    x  = circle[0]
    y= circle[1]
    r= circle[2]
    #if x - r >= -edges[0]/2 and y - r >= -edges[1]/2 and x+r <= edges[0]/2 and y+r<= edges[1]/2:
    if r <= x <= edges[0]-r and r <= y <= edges[1]-r:
        return True
    return False

weird_angles = []


@jit(nopython = True)
def get_disperse_end(circle ,edges , distance_faces , degree=5 , max_itera =200 , distribution = "normal"):

    
    if distribution == "normal":
        mean,var = degree , .5
        sampler=1
        # print(args.dispersion)
    elif distribution == "lognormal":
        mean ,var = WF.solve_for_log_normal_parameters(degree , 1)
        sampler = 0


    p = np.array([distance_faces , 0 , 0] )

    new_p = np.array( [edges[0]*2 , edges[0]*2 ,edges[0]*2 ] )
    val = is_inside_box(new_p ,edges)
    while max_itera > -1 and  not is_inside_box(new_p ,edges):
        max_itera -=1
        new_p   = WF.wattson_distribution(p , mean_dispersion = mean , var = var  , sampler = sampler)
        scale = distance_faces/new_p[0]
        new_p *= scale
        new_p = np.array([new_p[1] + circle[0], new_p[2] + circle[1] , circle[2]] )

    #weird_angles.append(angle)
    return new_p #, angle 




####3####3

def create_grid(L1 , L2 , l_split):
    n_ceil1 = np.ceil(L1/l_split).astype(int)
    n_ceil2 = np.ceil(L2/l_split).astype(int)
    
    grid_circles = [ [[] for _ in range(n_ceil2)] for _ in range(n_ceil1) ]
    return grid_circles , n_ceil1 , n_ceil2


def norml2(x,y):
    return np.sqrt(x**2 + y**2)

def norml2_single(x1,y1 , x2,y2):
    return np.sqrt((x1-x2)**2 + (y1-y2)**2)


def get_xy_minmax(x,y,r , l_split , n_ceil1, n_ceil2):

    x_min = np.floor((x-r)/l_split).astype(int)
    x_max = np.floor((x+r)/l_split).astype(int)
    y_min = np.floor((y-r)/l_split).astype(int)
    y_max = np.floor((y+r)/l_split).astype(int)
    diag = np.sqrt(2)*l_split
    
    x_min = np.max([0 , x_min])
    x_max = np.min([n_ceil1-1 , x_max])
    y_min = np.max([0 , y_min])
    y_max = np.min([n_ceil2-1 , y_max])

    return x_min , x_max , y_min , y_max , diag


def insert_circle_grid(grid_circles , new_circle , l_split , n_ceil1 , n_ceil2):
    x,y,r = new_circle

    x_min , x_max , y_min , y_max , diag = get_xy_minmax(x,y,r , l_split , n_ceil1 , n_ceil2)

    for xi in range(x_min , x_max+1):
        for yi in range(y_min , y_max+1):
            if norml2_single(xi*l_split , yi*l_split , x , y) < r+ 2*diag/2:
                grid_circles[xi][yi].append(new_circle)
    
    return

def search_intersection(grid_circles , new_circle , l_split , n_ceil1 , n_ceil2):
    x,y,r = new_circle

    xi = np.floor(x/l_split).astype(int)
    yi = np.floor(y/l_split).astype(int)

    x_min , x_max , y_min , y_max , diag = get_xy_minmax(x,y,r , l_split , n_ceil1, n_ceil2)

    for xi in range(x_min , x_max+1):
        for yi in range(y_min , y_max+1):
            for circle_i in grid_circles[xi][yi]:
                if norml2_single(circle_i[0] , circle_i[1] , x , y) < r+circle_i[2]:
                    return 1
    return 0





def circle_grid_magnitude(grid_circles):
    n = len(grid_circles)
    m = len(grid_circles[0])
    num_circles = np.zeros((n,m))
    for i in range(n):
        for j in range(m):
            num_circles[i,j] = len(grid_circles[i][j])
    #flip matrix so that when plotted, it is in the correct orientation
    return num_circles.T
###########



### given some points finds one point that doesn't intersect
def wrapper_tree_points22(points , grid_circle_middle , l_split , n_ceil1 , n_ceil2):
    """ tree_middle = pp[0] """
    """ dist , aa  = tree_middle.kneighbors( pp[1], 1, return_distance=True) """

    intersections = [ search_intersection(grid_circle_middle , point , l_split , n_ceil1 , n_ceil2) for point in points]

    idi_max = np.argmin(intersections)
    #print("suma no intersections" )
    return intersections[idi_max] , points[idi_max]


def split_pool_tree(points ,grid_circle_middle , l_split , n_cores , n_ceil1 , n_ceil2 , pool_obj = None):

    dists = []
    sub_points = []

    """ if len(points) > 40: """
    if len(points) > np.inf:
    #if len(points) > 10:
        splits = n_cores +1
        idxs = np.linspace(0,len(points), splits).astype(int)

        """ inputs = [ [ tree_middle , points[idxs[i]:idxs[i+1] ] ] for i in range(splits-1)]  """
        inputs = [   [points[idxs[i]:idxs[i+1]] , grid_circle_middle , l_split , n_ceil1 , n_ceil2  ] for i in range(splits-1)] 

        #print(len( inputs) , len(inputs[-1][1]))

        ##return tree_middle.kneighbors( inputs[0], 1, return_distance=True)
        
        start = time.time()
        for result in pool_obj.imap_unordered(wrapper_tree_points22, inputs):

            dists.append(result[0])
            sub_points.append(result[1])


        end = time.time()
    else:
        result = wrapper_tree_points22( points , grid_circle_middle , l_split , n_ceil1 , n_ceil2)

        dists.append(result[0])
        sub_points.append(result[1])
    #print(end - start , "'s elapsed")
    return dists , sub_points




counter_up =[]
counter_down =[]
angles_voxel = []
def get_gamma_packing_both_faces(rads, edges , distance_faces ,degree = 5, n_cores= 1 , l_split = .3):

    pool_obj = multiprocessing.Pool(n_cores)


    r = rads[0]
    ex ,ey = edges
    circle_down =np.array([np.random.uniform(r,ex-r), np.random.uniform( r,ey-r),r])
    #circle_down =np.array([np.random.uniform(-ex/2 + r,ex/2-r), np.random.uniform(-ey/2 + r,ey/2-r),r])
    circles_down = [circle_down]




    next_circle = get_disperse_end(circle_down , edges = edges , distance_faces= distance_faces , degree=degree)
    circles_up = [ next_circle ]

    grid_down , _ , _ = create_grid(ex , ey  , l_split)
    grid_up , n_ceil1 , n_ceil2   = create_grid(ex , ey , l_split)

    insert_circle_grid(grid_down , circle_down , l_split , n_ceil1 , n_ceil2)
    insert_circle_grid(grid_up   , next_circle , l_split , n_ceil1 , n_ceil2)
    
    #distas = []
    bar = progressbar.ProgressBar(max_value=len(rads) -1)
    for i , r in enumerate(rads[1:]) :
        if i%200 ==0:
            print("r = " , r )
        bar.update(i)
         
        """ tree_down = NearestNeighbors(n_neighbors=1 ,   n_jobs =1 , metric = my_distance , algorithm = "ball_tree") """
        """ tree_up   = NearestNeighbors(n_neighbors=1 ,   n_jobs =1 , metric = my_distance , algorithm = "ball_tree") """
        """"""
        """ tree_down.fit(circles_down) """
        """ tree_up.fit(circles_up) """


        settled_face = 0
        max_itera_down = 45 ## original number
        """ max_itera_down = 100 """
        while settled_face !=2 and  max_itera_down >0:
            max_itera_down -=1
            
            counter =10
            counter =5
            while 1:
                counter +=1
                new_points =[np.array([np.random.uniform( r,ex -r) , np.random.uniform(r,ey-r), r]) for a in range (10 + 10*counter)]
                #intersections, sub_points , = split_pool_tree( new_points ,n_cores ,tree_down , pool_obj)
                intersections,sub_points =  split_pool_tree(new_points , grid_down , l_split , n_cores , n_ceil1 , n_ceil2 , pool_obj = pool_obj )

                mini_id = np.argmin(intersections)
                if intersections[mini_id] == 0  :
                    settled_face =1 
                    break

            if settled_face ==0 :
                continue

            possible_down = sub_points[mini_id]

            max_itera_up = 100
            max_itera_up = 20
            counter =10
            counter =5
            while max_itera_up > 0:
                max_itera_up -=1
                counter +=1
                new_points =[get_disperse_end(possible_down , edges = edges , distance_faces= distance_faces , degree=degree)  for a in range (10 + 10*counter)]
                #dist, ind = tree_up.query(new_points, k=1) 

                #dist, ind , = tree_up.kneighbors( new_points, 1, return_distance=True)
                #dist, sub_points , = split_pool_tree( new_points ,n_cores, tree_up, pool_obj)

                intersections,sub_points =  split_pool_tree(new_points , grid_up , l_split , n_cores , n_ceil1 , n_ceil2 , pool_obj = pool_obj)
                #max_zero = np.where(dist >=max_limit)[0]
                #if len(max_zero) ==0:
                #    continue
                #mini_id = max_zero[dist[max_zero].argmin()]

                mini_id = np.argmin(intersections)
                if intersections[mini_id] == 0  :
                    settled_face =2

                    #counter_up.append(counter)
                    break

            if settled_face ==2:
                possible_up = sub_points[mini_id]
                circles_down.append(possible_down)
                circles_up.append(possible_up)
                insert_circle_grid(grid_down , possible_down , l_split , n_ceil1 , n_ceil2)
                insert_circle_grid(grid_up   , possible_up , l_split , n_ceil1 , n_ceil2)
            else: 
                settled_face=0


        if max_itera_down == 0 and settled_face !=2:
            ## if this, error while packing the circles
            print ("didn't fit circles TOO PACKED!")


    pool_obj.close()
    mag1 = circle_grid_magnitude(grid_down )

    mag2 = circle_grid_magnitude(grid_up )

    #print("magnitudes grids: ")
    """ plt.subplot(1,2,1) """
    """ plt.imshow(mag1) """
    """ plt.subplot(1,2,2) """
    """ plt.imshow(mag2) """

    #center center circles in square -ex/2 , ex/2 , -ey/2 , ey/2

    for i in range(len(circles_down)):
        circles_down[i][0] -= ex/2
        circles_down[i][1] -= ey/2
        circles_up[i][0] -= ex/2
        circles_up[i][1] -= ey/2

    return circles_down , circles_up


def into_3dspace(positions , axis  , val , offset= [0,0,0]):
    new_pos = []
    for xi,yi , ri in positions:
        c2 = [xi,yi]
        i2 = 0
        i3 = 0
        c3 = np.array([0.0]*3)
        while i3 < 3:
            if i3== axis:
                c3[i3] =val + offset[i3]
            else:
                c3[i3] = c2[i2] + offset[i3]
                i2+=1
            i3+=1
        new_pos.append([*c3 ,ri])
    new_pos = np.array(new_pos) 
    return new_pos

  
def Rx(theta):
  return np.matrix([[ 1, 0           , 0           ],
                   [ 0, m.cos(theta),-m.sin(theta)],
                   [ 0, m.sin(theta), m.cos(theta)]])
  
def Ry(theta):
  return np.matrix([[ m.cos(theta), 0, m.sin(theta)],
                   [ 0           , 1, 0           ],
                   [-m.sin(theta), 0, m.cos(theta)]])
  
def Rz(theta):
    return np.matrix([[ m.cos(theta), -m.sin(theta), 0 ],
                   [ m.sin(theta), m.cos(theta) , 0 ],
                   [ 0           , 0            , 1 ]])

def alling_bundle_with_vector(face1 , face2 ,  theta = 45):

    #my_mat  =rotation_matrix_from_vectors(axis , new_axis)
    my_mat  = Rz(WF.degree_to_radians(theta))

    rad1 = face1[:,-1]
    rad2 = face2[:,-1]

    mid_point = np.average(face1)

    f1= my_mat.dot(face1[:3, :])
    f2= my_mat.dot(face2[:3, :])

    face1[:3 , :] = f1
    face2[:3 , :] = f2

    return face1.T , face2.T
 

def analysis_dispersion_axons_start_end(bundle1_start , bundle1_end , vx):
    angles = []
    for c1 , c2 in zip(bundle1_start , bundle1_end):
        vector = np.array(c1[0:3]) - np.array(c2[0:3])
        my_angle = get_angle(vector , vx)
        angles.append(my_angle)


    #plot a gaussian density function with mean mu and standard deviation sigma
    mu, sigma = args.dispersion , .5 # mean and standard deviation
    print("variables" , mu , sigma)
    x = np.linspace(-1, 1, 100)

    #sample n points from gaussiang with mu and sigma
    n = 1000
    samples = np.random.normal(mu, sigma, n)

    #plot the histogram as density of samples
    #sns.histplot(samples, stat='density')

    sns.histplot(angles , stat='density' , label = 'Empirical density')

    sns.kdeplot(samples , bw_adjust=2 , label = 'True density')
    sns.kdeplot(angles , bw_adjust =2 ,   label = 'Empirical density')

    plt.axvline(x = np.mean(angles) , c = 'r' , label = 'Empirical mean angle')
    plt.legend()

    print("mean angle , "  , np.mean(angles))
    if args.verbose:
        plt.show()


def get_radii_distribution_from_rect_icvf(theta = 1.1 , kappa = .5 , edges = [10 ,10] , icvf = .5  , rad_threshold = [.5,4]  , bias = 0 ):
    area_sum = 0
    max_area =  edges[0]*edges[1]*icvf

    rads_b1= [] 
    while  area_sum < max_area:
        r1 = np.random.gamma(theta, kappa,1)[0]
        #r1 -= bias
        if rad_threshold[1] >=  r1 >= rad_threshold[0]:
            r1 += bias
            rads_b1.append(r1)
            area_sum += np.pi*r1*r1

    rads_b1 = sorted(rads_b1)[::-1] 
    return np.array(rads_b1)


#%%


def create_voxel(args):


    if type(args.gamma_theta) == list:
        print("gamma_pairs found")
        args.gamma_theta , args.gamma_kappa = args.gamma_pairs

    # print long line of dashes
    print("Starting voxel generation voxel")
    print("-"*80)
    print(args)
    print("************")

    np.random.seed(args.number_of_repetitions)


    small_buffer = 0
    lenght1 = args.lenght_side + small_buffer ### plus small buffer

    if args.n_bundles ==2:
        crossing_angle = WF.degree_to_radians(args.crossing_angle)
        overlap = args.overlap
        intersection_lenght = overlap*lenght1
        small_edge = (lenght1 + intersection_lenght)/2

    elif args.n_bundles ==1:
        small_edge = lenght1
    else: 
        raise Exception("Wrong number of bundles")




    icvf = args.icvf






    edges1 = np.array([float(lenght1), float(small_edge)])
    ic(edges1)
    # edges2 =np.array( [-1.0 , -1.0])

    rad_threshold =  [args.min_rad + args.bias , args.max_rad + args.bias]
    rads_b1 = get_radii_distribution_from_rect_icvf( theta = args.gamma_theta , kappa = args.gamma_kappa, edges = edges1 , rad_threshold = rad_threshold , icvf = icvf , bias = args.bias)
    mean_radb1 = np.mean(rads_b1)
    """ rads_b1 = get_radii_distribution_from_rect_icvf( theta = 2.8 , kappa = .13, edges = edges1 , rad_threshold = [.25,2] , icvf = icvf , bias = args.bias -.25) """
    




    if args.n_bundles==2:
        # lenght2 = np.sin(crossing_angle)*( lenght1 + lenght1/np.tan(crossing_angle))
        lenght2 = ( np.sqrt(2) -1 )*np.sin(2*crossing_angle/180*np.pi)   *lenght1 +  lenght1

        edges2 = [float(lenght2), float(small_edge)]
        ic(edges2)
        rads_b2 = get_radii_distribution_from_rect_icvf( theta = args.gamma_theta , kappa = args.gamma_kappa, edges = edges2 , rad_threshold = rad_threshold , icvf = icvf)

        print(f"L = {edges1[0]} x {edges1[1]} ,    L2 = {edges2[0]} x {edges2[1]}")

#plt.subplot(1,2,1)
    if args.verbose:
        """ plt.hist(rads_b1 - args.bias , color = 'b' , density = True , label = "true dist" , alpha = 0.5 , bins =100) """
        plt.hist(rads_b1  , color = 'r' , density = True , label = "biassed dist" , alpha = 0.5 , bins =50)
        plt.axvline(x = np.mean(rads_b1) , c = 'k' , label = "mean rad {0:.4f}".format(np.mean(rads_b1)))
        """ plt.axvline(x = np.mean(rads_b1) + args.bias , c = 'k' , linestyle = '--' , label = "mean rad {0:.4f}".format(np.mean(rads_b1) + args.bias)) """
        plt.legend()
        #xrange
        plt.xlim([-0.05, np.max(rads_b1)])
        ic(edges1)
        rad0 =  rads_b1-args.bias
        #sort from biggest to smallest
        rad0 = sorted(rad0)[::-1]
        radii_file = args.outfile.split(".")[0] + "_radii.txt"
        #radii_file2 = args.outfile2.split(".")[0] + "_radii.txt"
        np.savetxt(radii_file , rad0)
        #np.savetxt(radii_file2 , rad0)



        if args.n_bundles==2:
            plt.hist(rads_b2 , color = 'r' , bins = 100)
            print(f"L = {edges1[0]}   L2 = {edges2[0]}")
        plt.title("Rad histogram histograms")

        plt.show()

#%%
    if args.n_bundles ==2:
        """ height1 =  lenght1 + 2*lenght1/np.tan(crossing_angle) """
        """ height2 =  lenght1/np.sin(crossing_angle) + lenght2/np.tan(crossing_angle) """

        height1 = args.depth_lenght_bundle
        height2 = args.depth_lenght_bundle

        ic(height1 ,height2)

    else:
        height1 = lenght1 
        height2 = -1 
    print(f" H1 = {height1}   H2 = {height2}")
    


#%%

    print("------------------------")
    print("Fitting starting and ending points ...  ")
    print(args.dispersion)
    print("------------------------")
    face_b1 , face_b11 = get_gamma_packing_both_faces(rads_b1, edges1 , height1 ,degree = args.dispersion , n_cores = args.n_cores , l_split = mean_radb1*4.5)


#plt.hist(weird_angles)
#plt.title("true angles")

    if args.n_bundles ==2:
        face_b2 , face_b22 = get_gamma_packing_both_faces(rads_b2, edges2 , height2 ,degree = args.dispersion , n_cores =args.n_cores , l_split = mean_radb1*1.5)
#print(len(rads_b1) , len(face_b1))






#%%
    if args.verbose :
        plt.subplot(args.n_bundles,2,1)
        plot_circle_array(face_b1)
        plt.subplot(args.n_bundles,2,2)
        plot_circle_array(face_b11)

        if args.n_bundles ==2:
            plt.subplot(2,2,3)
            plot_circle_array(face_b2)
            plt.subplot(2,2,4)
            plot_circle_array(face_b22)
        plt.show()

#%%

## define one population    ## seeting in coordinates z -> 2
    if args.n_bundles ==1:
        face3d_b1 = into_3dspace(face_b1 , 0    ,-height1/2 , offset = [0,0,0])
        face3d_b11   = into_3dspace(face_b11, 0 , height1/2 , [0 , 0,0])

        """ face3d_b1 = into_3dspace(face_b1 , 0    , 0, offset = [0,0,0]) """
        """ face3d_b11   = into_3dspace(face_b11, 0 , height1 , [0 , 0,0]) """

    elif args.n_bundles ==2:

        offset_intersecteion = (lenght1 - small_edge)/2
        face3d_b1 = into_3dspace(face_b1    , 0 ,-height1/2 , offset = [0 , 0,  offset_intersecteion])
        face3d_b11   = into_3dspace(face_b11, 0 , height1/2 , offset = [0 , 0,  offset_intersecteion])

        face3d_b2 = into_3dspace(face_b2    , 0 ,-height2/2 , offset = [0 , 0, -offset_intersecteion])
        face3d_b22   = into_3dspace(face_b22, 0 , height2/2 , offset = [0 , 0, -offset_intersecteion])
    else :
        raise Exception("Wrong number of bundles")
#%%

#new_axis1 = np.repeat(1,3)
    new_axis1 = np.repeat(0,3)
    new_axis1[0] = 1 




    if args.n_bundles ==2:
        theta = crossing_angle
        phi = WF.degree_to_radians(0)  # np.random.uniform(0  , 2*np.pi)
        new_axis2 = WF.polar_to_cartesian(phi ,  theta , 1)
        new_axis2 = new_axis2/np.linalg.norm(new_axis2)



    face3d_b1 ,face3d_b11 = alling_bundle_with_vector(face3d_b1.T ,face3d_b11.T , theta = 0)

    z = np.array([1,0,0])
    all_angles = []
    
    if args.verbose:
        for p1,p2 in zip(face3d_b1 , face3d_b11):
            v = p2 - p1
            v=v[:3]
            v = v/np.linalg.norm(v)
            angles = np.arccos(np.dot(v,z))*180/np.pi
            all_angles.append(angles)
        plt.hist(all_angles)
        plt.axvline(np.mean(all_angles) , color = 'r')
        plt.title("angles")
        plt.show()





    if args.n_bundles ==2:
        face3d_b2 ,face3d_b22 = alling_bundle_with_vector(face3d_b2.T ,face3d_b22.T , theta = args.crossing_angle)


#%%

    # analysis_dispersion_axons_start_end(face3d_b1 , face3d_b11  ,new_axis1)
    # if args.n_bundles ==2:
    #     analysis_dispersion_axons_start_end(face3d_b2 , face3d_b22  ,new_axis2)

#%%

    if args.n_bundles ==2:
        starts = np.concatenate([face3d_b1 , face3d_b2])
        ends = np.concatenate([face3d_b11 , face3d_b22])
    else: 
        starts = face3d_b1
        ends = face3d_b11

    cilindros = []
    cilindros_sinus = []

    n_controls = []
    print("creating streamlines")
    bar = progressbar.progressbar
    for si , ei in bar(zip( starts , ends)):
        ##print(si , ei)

        """ if args.wavy_radii: """
        """     cil  , n_control = RWS.create_streamline(si , ei , PerExtraDist=-.1 , wavy_radii = args.wavy_radii  , lenght_stick=5 , lenght_ball= 5 , min_rad = args.min_rad , max_rad = args.max_rad) """
        """     cilindros_sinus.append(cil) """


        cil  , n_control = RWS.create_streamline(si , ei , PerExtraDist=-.1 , wavy_radii = False)
        n_controls.append(n_control)

        cilindros.append(cil)
#%%

    if args.verbose:

        plt.hist(n_controls)
        plt.title("Histogram Control points in the strands")
        plt.show()

#%%

#write_jony(cilindros , "start_1bundle.txt")
    print(f"... Saving in {args.outfile}")
    RWS.write_strand_list_with_box(cilindros       , args.outfile , box_size = args.lenght_side)
    """ RWS.write_strand_list_with_box_jony(cilindros       ,  "outer_" + args.outfile , box_size = args.lenght_side , scale_rad =1) """
    """ RWS.write_strand_list_with_box_jony(cilindros       ,  "inner_" + args.outfile , box_size = args.lenght_side ,scale_rad =.7) """

    """ if args.wavy_radii: """
    """     RWS.write_strand_list_with_box(cilindros_sinus , args.outfile2, box_size = args.lenght_side) """

    
    RC.write_my_config(args)
    ## write file configuratoins in the same place as the outfile





if __name__ == "__main__":

    if args.config_file is not None:
        args = RC.read_config_file(args.config_file)


        

    print(" ------------------ ")
    print(args)
    print(" ------------------ ")

    np.random.seed()


    if args.n_bundles>2:
        print("three bundles not implemented yet")
        sys.exit()
    create_voxel(args)




