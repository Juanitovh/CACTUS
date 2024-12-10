

###
import numpy as np 
import matplotlib.pyplot as plt
from rich import print
from rich.progress import track
from rich.console import Console
from icecream import ic
###

import scipy
#len_strand = cylinder[0][0:3] - cylinder[-1][0:3]


def fit_fat_ball( lenght_ball , plot=False):
    #modifucation of sine wave, to be fatter
    f2 = .5
    pos1d = np.arange(0, 1.0001, .01)
    pow = 3.5
    scale_tail = 0.65
    sinewave1 = np.abs(np.sin(pos1d * 2 * np.pi *f2 + np.pi/4 + np.pi/16 ) )**(pow+1)
    sinewave1= sinewave1 / np.max(sinewave1)*scale_tail
    sinewave2 = np.abs(np.sin(pos1d * 2 * np.pi *f2 ) )  **pow
    sinewave2= sinewave2 / np.max(sinewave2)
    sinewave3 = np.abs( np.sin(pos1d * 2 * np.pi *f2 + np.pi/4*3 - np.pi/16 )  ) **(pow+1)
    sinewave3= sinewave3 / np.max(sinewave3)*scale_tail

    sum_sinewaves = sinewave1 + sinewave2 + sinewave3
    sum_sinewaves = sum_sinewaves / np.max(sum_sinewaves)

    minimum = np.min(sum_sinewaves)
    sum_sinewaves = sum_sinewaves - minimum
    maximum = np.max(sum_sinewaves)
    sum_sinewaves = sum_sinewaves / maximum



    pos_real_ball = pos1d * lenght_ball
    y_interp = scipy.interpolate.interp1d(pos_real_ball , sum_sinewaves, kind='cubic')
    sum_sinewaves_interpol = y_interp(pos_real_ball)

    if plot :

        plt.plot(pos1d, sinewave2**1)
        plt.plot(pos1d, sinewave2**3)
        plt.plot(pos1d, sinewave2**4)
        plt.plot(pos1d, sinewave2**10)

        plt.show()
        plt.plot(pos1d, sinewave1)
        plt.plot(pos1d, sinewave2)
        plt.plot(pos1d, sinewave3)
        plt.plot(pos1d, sum_sinewaves ,"k--")
        plt.plot(pos1d, sum_sinewaves_interpol+.01 ,"r--")

    return y_interp


import sys
def oscillating_radii( lenght_strand, n_control_points, lenght_ball, lenght_stick, plot=False):

    offset = np.random.uniform(0 , lenght_ball + lenght_stick)
    #offset = lenght_ball
    offset=11

    y_interp = fit_fat_ball( lenght_ball=lenght_ball, plot=plot)


    id_capsule = np.linspace(offset, lenght_strand+offset , n_control_points)
    step = id_capsule[1] - id_capsule[0]
    final_radii = [] 
    for p in id_capsule:
        p = p% (lenght_ball + lenght_stick )
        if  0 <=p and  p < (lenght_ball + 0):
            final_radii.append(y_interp(p ) +1)
        else:
            final_radii.append(1)

    n_steps_ball = int(lenght_ball / step)
    final_radii = np.array(final_radii)

    if np.mean(final_radii[:n_steps_ball]) !=1:
        final_radii[:n_steps_ball] = 1
        
        ii = n_steps_ball
        while final_radii[ii] != 1:
            final_radii[ii] = 1
            ii += 1
    if np.mean(final_radii[-n_steps_ball:]) !=1:
        final_radii[-n_steps_ball:] = 1
        ii = -n_steps_ball-1
        while final_radii[ii] != 1:
            final_radii[ii] = 1
            ii -= 1





    final_radii = np.array(final_radii)
    return id_capsule - offset, final_radii

if __name__ == "__main__":
    id_capsule, final_radii = oscillating_radii( lenght_strand=35, n_control_points = 500, lenght_ball=5, lenght_stick=5, plot=True)

    plt.show()
    plt.plot(id_capsule, final_radii)
    plt.show()












