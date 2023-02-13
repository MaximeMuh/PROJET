# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 15:27:40 2023

@author: maxim
"""
import math
import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate
#paramètres pour siumuler des lois normales
nu=10
sigma=4
#parametres pour la simulation et la chaîne de Markov
dt=0.1
#on conseille de prendre N>=10^5


#on définit plusieurs fonctions tests que l'on veut échantilloner
def f(x):
    return np.exp(-((x-nu)**2)/(2*(sigma**2)))
def f2(x):
    return np.exp(-(abs(x)))
def V(x):
    return np.exp(-np.cos(2*np.pi*x))

#on définit notre fonction instrumentale, on prend ici une gaussienne
def q(x,x_n):
    return (1/math.sqrt(2*math.pi))*math.exp(-((x-x_n)**2)/2)


#on définit la fonction acceptation qui accepte ou rejette une proposition selon sa valeur
def acceptation(prop,x_n,fonct):
    alpha=min(1,(fonct(prop)*q(x_n,prop))/(fonct(x_n)*q(prop,x_n)))
    u=np.random.uniform(0,1)
    if u<=alpha:
        return prop
    else:
        return x_n



#on définit la foncion suivant la méthode Metropolis-Hastings ou l'on calcule les élements de notre suite de Markov, on trace ensuite un histogramme que l'on compare avec la fonction test
def MH(x0,N,fonct):
    l=[x0]
    for k in range (N-1):
        x0=np.random.normal(0,1)*np.sqrt(2*dt)+l[-1]
        l.append(acceptation(x0,l[-1],fonct))
    l=np.array(l)
    if (fonct==V):
        l=abs(l)%1
        plt.hist(l,bins=50,density=True ,color="lightblue",edgecolor='blue')
        plt.title("Echantillonage selon MH pour V(x)=cos(2*π*x),N=10^6,dt=0.1")
        l2=np.arange(0,1,0.001)
        A=scipy.integrate.quad(fonct,0,1)[0]
    else:
        if fonct==f:
            l2=np.arange(-5*sigma+nu,nu+5*sigma,0.001)
            plt.title("Echantillonage selon MH d'une gausienne,nu=10,sigma=4,N=10^6,dt=0.1")
        else:
            l2=np.arange(-5,5,0.001)
            plt.title("Echantillonage selon MH , N=10^6,dt=0.1")
        plt.hist(l,bins=50,density=True ,color="lightblue",edgecolor='blue')
        A=scipy.integrate.quad(fonct,-100,100)[0]
    y=fonct(l2)/A
    plt.plot(l2,y,color='green')
    plt.show()




