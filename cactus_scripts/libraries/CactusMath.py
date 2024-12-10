###
import numpy as np 
import matplotlib.pyplot as plt
from rich import print
from rich.progress import track
from rich.console import Console
from icecream import ic
###


def exponential_linspace(start , end , num , exp=3):
    line = ( np.linspace(0,1,num , endpoint=True))**exp
    line =  ( line - line[0])/ (line[-1]- line[0])
    
    line = end- (end - start)*line+start
    return line[::-1]

def exponential_linspace_int(start , end , num , exp=3):
    ans = exponential_linspace(start , end , num , exp=exp)
    ans =ans.astype(int)
    return np.unique( ans)



