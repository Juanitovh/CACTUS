###
import numpy as np 
import matplotlib.pyplot as plt
from rich import print
from rich.progress import track
from rich.console import Console
from icecream import ic
###
#clear temporal files python clear ram dump
import os
import gc
import time
import shutil


def clear_ram_dump():
    gc.collect()
    print("Ram Dumped")


def exponential_linspace(start , end , num , exp=3):
    line = ( np.linspace(0,1,num , endpoint=True))**exp
    line =  ( line - line[0])/ (line[-1]- line[0])
    
    line = end- (end - start)*line+start
    return line[::-1]

def exponential_linspace_int(start , end , num , exp=3):
    ans = exponential_linspace(start , end , num , exp=exp)
    ans =ans.astype(int)
    return np.unique( ans)



def paste_half(x):
    ## x must be a list of n elements containint two-ples of (local itera dim)
    n = len(x)
    x2 = []
    print("Binary merging of ", n, " elements")
    for i in track(range(n//2)):
        x2.append( x[2*i] + x[2*i+1]  )
    if n%2 ==1:
        x2.append(  x[-1])
    del x
    return x2

def paste_binary(x):
    while len(x)>=2:
        x = paste_half(x)
        gc.collect()
    return x[0]
