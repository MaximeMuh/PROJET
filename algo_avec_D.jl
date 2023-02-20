# On importe les bibliotheques
using Random, Distributions, Plots, SciPy

# On definit les constantes et les fonctions utilisees dans la suite
N = Int32(1e6) # Nombre de configurations calculees
dt = 1e-3 # Pas de temps
frequence = 1 # Frequence du potentiel V (qui est un cosinus)
nb_bins = 100 # Nombre de bins pour l'affichage de l'histogramme 
              # et nombre de point pour le trace de mu

V(x) = cos(2 * pi * frequence * x) # Potentiel
grad_V(x) = -2 * pi * sin(2 * pi * frequence * x)
D(x) = 1.5 + cos(2 * pi * x)
div_D(x) = -2 * pi * sin(2 * pi * x)
U = Uniform() # Loi uniforme
G = Normal() # Gaussienne centree reduite

# attention : 
log_q(x, y) = - (x - y + dt * (D(y) * grad_V(y) - div_D(y)))^2 / (4 * D(y) * dt) - log(D(y)) / 2

# Lq contiendra la suite finie des configurations
Lq = zeros(Float32, N)

for k = 1:(N-1)
    q_k = Lq[k]

    # On calcule la nouvelle configuration a partir de l'equation de Langevin
    candidat = q_k + (-D(q_k) * grad_V(q_k) + div_D(q_k))* dt + sqrt(2 * D(q_k) * dt) * rand(G)

    # Metropolis-Hastings, on calcule log_alpha car c'est plus rapide   
    log_alpha = min(0, V(q_k) - V(candidat) + log_q(q_k, candidat) - log_q(candidat, q_k))
                    
    if log(rand(U)) <= log_alpha
        Lq[k + 1] = candidat
    else
        Lq[k + 1] = q_k
    end
end

# On ramene toutes les valeurs de Lq dans le tore
for k = 1:N
    Lq[k] = abs(Lq[k]) % 1
end


## Affichage

# On definit la densite mu de la loi qu'on veut simuler, pour cela on va
# d'abords calculer le facteur de normalisation A
A = SciPy.integrate.quad(x -> exp(-V(x)), 0, 1)[1]
mu(x) = exp(-V(x)) / A

# On realise le graphique et on l'affiche
x = range(0, 1, length = nb_bins)
histogram(Lq, bins = nb_bins, normalize=:true, label = "simulation dt=$dt N=$N",
title = "Echantillonage pour V(x)=cos($(2 * frequence)*pi*x)")
display(plot!(x, mu, label = "mu"))