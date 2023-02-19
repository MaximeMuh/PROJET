# On importe les bibliotheques
using Random, Distributions, Plots, QuadGK

# On definit les constantes et les fonctions utilisees dans la suite
N = Int32(1e6) # Nombre de configurations calculees
dt = 1e-4 # Pas de temps
frequence = 1 # Frequence du potentiel V (qui est un cosinus)
nb_bins = 100 # Nombre de bins pour l'affichage de l'histogramme 
              # et nombre de point pour le trace de mu


function D(x)
    exp(x)
end 
function q(x,y)
    exp((-1)/(4*dt)*(x-y-dt*deriv_V(x))^2)
end 
V(x) = cos(2 * pi * frequence * x)
deriv_V(x) = -2 * pi *frequence*  sin(2 * pi * frequence * x)
deriv_D(x)= exp(x)
U = Uniform()
G = Normal()

# Lq contiendra la suite finie des configurations
Lq = zeros(Float32, N)
Lq[1] = rand()

for k = 1:(N-1)
    # On calcule la nouvelle configuration a partir de l'equation de Langevin
    candidat = Lq[k] + (-D(Lq[k]) * deriv_V(Lq[k]) + deriv_D(Lq[k])) * dt + sqrt(2 * D(Lq[k]) * dt) * rand(G)


    # Metropolis-Hastings
    alpha= min(1,(exp(-V(candidat))*q(candidat, Lq[k])/ exp(-V(Lq[k]))*q(Lq[k], candidat)))

    if (rand(U)) <= alpha
        Lq[k + 1] = candidat
    else
        Lq[k + 1] = Lq[k]
    end
end

# On ramene toutes les valeurs de Lq dans le tore
    for k= 1:N
        Lq[k]= Lq[k]%1
end 


## Affichage

# On definit la densite mu de la loi qu'on veut simuler, pour cela on va
# d'abords calculer le facteur de normalisation A
A = quadgk(x -> exp(-V(x)), 0, 1)[1]
mu(x) = exp(-V(x)) / A

# On realise le graphique et on l'affiche
x = range(0, 1, length = nb_bins)
histogram(Lq, bins = nb_bins, normalize=:true, label = "simulation",
title = "Echantillonage pour V(x) = cos($(2 * frequence)*pi*x)")
display(plot!(x, mu, label = "mu"))


