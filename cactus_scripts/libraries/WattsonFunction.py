import numpy as np
from numba import jit

@jit(nopython = True)
def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    #if distance between vectors is zero, return identity matrix
    if np.linalg.norm(vec1 - vec2) < .001:
        return np.identity(3)
    
    a, b = (vec1 / ((np.linalg.norm(vec1) ) )).reshape(3), (vec2 / ((np.linalg.norm(vec2) ))) .reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

@jit(nopython = True)
def degree_to_radians(angle):
    return angle*np.pi/180

@jit(nopython = True)
def radians_to_degree(angle):
    return angle*180/np.pi

@jit(nopython = True)
def polar_to_cartesian(phi,theta , r):
    p = np.zeros(3)
    p[0] = r* np.cos(phi)*np.sin(theta)
    p[1] = r* np.sin(phi)*np.sin(theta)
    p[2] = r*np.cos(theta)
    return p

@jit(nopython = True)
def my_norm(p):
    suma = 0
    for pi in p:
        suma += pi*pi
    return np.sqrt(suma)

@jit(nopython = True)
def wattson_distribution_bad(my_vec , mean_dispersion =5 , distribution = 0):
    vecz = np.zeros(3)
    vecz[2] = 1

    my_vec = my_vec/my_norm(my_vec)

    m_rot = rotation_matrix_from_vectors(vecz , my_vec)

    #std = np.sqrt(np.pi/2)*mean_dispersion
    #angle = np.abs(np.random.normal(0 , std))
    angle = np.abs( np.random.normal(mean_dispersion, .5))
    #angle = np.abs(np.random.normal(0 , mean_dispersion))
    theta = degree_to_radians(angle)
#    print(theta)
#    theta = np.arccos(theta)
    phi = np.random.uniform(0,2*np.pi)
    p = polar_to_cartesian(phi ,theta , 1)

    p = m_rot.dot(p.T).T
    return p 

from numpy import log

@jit(nopython= True, fastmath = True)
def solve_for_log_normal_parameters(mean, variance):
    if mean ==0:
        return 0,1
    sigma2 = log(variance/mean**2 + 1)
    mu = log(mean) - sigma2/2
    return (mu, sigma2)

@jit(nopython = True)
def wattson_distribution(my_vec , mean_dispersion =5, var=1 ,sampler = 0):
    vecz = np.zeros(3)
    vecz[2] = 1

    my_vec = my_vec/my_norm(my_vec)

    m_rot = rotation_matrix_from_vectors(vecz , my_vec)

    #std = np.sqrt(np.pi/2)*mean_dispersion
    #angle = np.abs(np.random.normal(0 , std))
#    if distribution =="normal":

    if mean_dispersion == 0 :
        angle = 0
    else :
        if sampler ==0:
            angle = np.random.lognormal(mean_dispersion, var)
        else :
            angle = np.abs(np.random.normal(mean_dispersion , var))



 #   elif distribution == "uniform":
 #       angle = np.random.uniform(0, mean_dispersion*2)
  #  elif distribution == "parallel":
   #     angle = 0
    #else :
     #   print("unavailable distribution")
      #  return
    #angle = np.abs(np.random.normal(0 , mean_dispersion))
    theta = degree_to_radians(angle)
#    print(theta)
#    theta = np.arccos(theta)
    phi = np.random.uniform(0,2*np.pi)
    p = polar_to_cartesian(phi ,theta , 1)

    p = m_rot.dot(p.T).T
    return p 
#%%
import matplotlib.pyplot as plt
#print("before __name__ guard")
if __name__ == '__main__':
    dispersion = 25

    theta_25 = degree_to_radians(25)
    v_25 =  np.array([0,np.tan(theta_25),1])
    v_25 = v_25/my_norm(v_25)

    my_vec = np.array([1,1,1])/np.sqrt(3)
    my_vec = np.array([0,0,1])
    distance_angles = []

    points = []
    angles = []
    for j in range(500):
        p  = wattson_distribution(my_vec , mean_dispersion = dispersion, var = 5 , sampler = 1)
        points.append(p)
        distance_angles.append(radians_to_degree(np.arccos(np.dot(v_25,p)/my_norm(p))))

        #angles.append(alfa)

    points = np.array(points).T

    plt.hist(distance_angles )
    mean = np.mean(distance_angles)
    plt.axvline(mean , color = "red")
    plt.title("angle distribution")
    plt.show()
    #print(np.mean(angles))
    #plt.hist(angles)
    #plt.show()

    import matplotlib.pyplot as plt 

    ax = plt.axes(projection='3d')

# Data for a three-dimensional line
    x,y,z = points
    ax.scatter3D(x, y, z);
    ax.scatter3D(0, 0, 0 , 'black' , s = 200);
    ax.scatter3D(*my_vec  , 'red' , s = 100 );
    ax.scatter3D(*v_25  , 'green' , s = 700 );
    plt.show()
