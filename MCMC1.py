# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 15:27:40 2023

@author: maxim
"""
import math
import numpy as np
from matplotlib import pyplot as plt
def f(x):
    return (1/np.sqrt(2*np.pi))*np.exp(-(x**2)/2)
def f2(x):
    return np.exp(-(abs(x)))/2


def q(x,x_n):
    return (1/math.sqrt(2*math.pi))*math.exp(-((x-x_n)**2)/2)



def acceptation(prop,x_n):
    alpha=min(1,(f(prop)*q(x_n,prop))/(f(x_n)*q(prop,x_n)))
    u=np.random.uniform(0,1)
    if u<=alpha:
        return prop
    else:
        return x_n




def MH(x0,N):

    l=[x0]
    for k in range (N-1):
        x0=np.random.normal(l[-1],1)
        l.append(acceptation(x0,l[-1]))

    l=np.array(l)
    print(l)
    plt.hist(l,bins=50,density=True ,color="lightblue")
    l2=np.arange(-7,7,0.001)
    y=f(l2)
    plt.plot(l2,y)
    plt.show()

