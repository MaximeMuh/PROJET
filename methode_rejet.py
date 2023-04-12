#On souhaite simuler une loi de probabilité de densité p. On suppose qu'on sait simuler une loi de proba de densité g tel que p <= k*g 
#où k est une constante
#On se donne sur l'univers une suite (Yi, Ui) de v.a IID avec Y1 de densité g et U1 uniforme sur 0:1
import random
from matplotlib import pyplot as plt
import numpy as np

def methode_rejet(p, g, N):
     echantillonx = []
     echantillony = []
     echantillonx_rejete = []
     echantillony_rejete = []
     while len(echantillonx) < N:
        # On prend un point compris entre 0 et 1 par la loi uniforme
        x = random.random()  
        #On calcule une hauteur aléatoire relative à la fonction qui borne p
        y = random.random() * g(x)     
        # On garde l'échantillon si la hauteur est en dessous de la courbe de p au même point
        if y <= p(x):
            echantillonx.append(x)
            echantillony.append(y)
        else:
            echantillonx_rejete.append(x)
            echantillony_rejete.append(y)
     echantillonx = np.array(echantillonx)
     echantillony = np.array(echantillony)
     echantillonx_rejete = np.array(echantillonx_rejete)
     echantillony_rejete = np.array(echantillony_rejete)
     return echantillonx, echantillony, echantillonx_rejete, echantillony_rejete
     
def affiche_simulation(p,g,N):
    x = np.linspace(0,1,5000)
    yp = p(x)
    yg = g(x)
    xsim, ysim, xrate, yrate = methode_rejet(p,g,N)
    plt.scatter(xsim, ysim, s=10, c='green', marker = '.', label='points gardés')
    plt.scatter(xrate, yrate, s=10, c='red', marker = '.', label='points rejetes')
    plt.plot(x, yp, "y-", linewidth= 3, label = 'densite p')
    plt.plot(x, yg, "m-" ,linewidth = 3, label = 'densite g')
    plt.show()

#TEST

def p(x):
    return  x*x
def g1(x):
    return x+1 - x
q = g1
N=500

affiche_simulation(p,g1,N)