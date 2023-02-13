# On importe les bibliotheques
using Random, Distributions, Plots

# On definit les constantes et les fonctions utilisees dans la suite
N = Int32(1e6) # Nombre de configurations calculees
dt = 1e-3 # Pas de temps
frequence = 1 # Frequence du potentiel V (qui est un cosinus)
nb_bins = 100 # Nombre de bins pour l'affichage de l'histogramme 
              # et nombre de point pour le trace de mu

V(x) = cos(2 * pi * frequence * x)
grad_V(x) = -2 * pi * sin(2 * pi * frequence * x)
U = Uniform()
G = Normal()

# Lq contiendra la suite finie des configurations
Lq = zeros(Float32, N)
Lq[1] = 0.5

for k = 1:(N-1)
    # On calcule la nouvelle configuration a partir de l'equation de Langevin
    candidat = Lq[k] - grad_V(Lq[k]) * dt + sqrt(2 * dt) * rand(G)

    # Metropolis-Hastings, on calcule log_alpha car c'est plus rapide
    log_alpha = min(0, V(Lq[k]) - V(candidat) + ((candidat - Lq[k] + dt * grad_V(Lq[k]))^2
    - (Lq[k] - candidat + dt * grad_V(candidat))^2) / (4 * dt))
    if log(rand(U)) <= log_alpha
        Lq[k + 1] = candidat
    else
        Lq[k + 1] = Lq[k]
    end
end

# On ramene toutes les valeurs de Lq dans le tore
for k = 1:N
    Lq[k] = abs(Lq[k]) % 1
end


## Affichage

# On definit la densite mu de la loi qu'on veut simuler, pour cela on va
# d'abords calculer le facteur de normalisation A
A = 1.266
mu(x) = exp(-V(x)) / A

# On realise le graphique et on l'affiche
x = range(0, 1, length = nb_bins)
histogram(Lq, bins = nb_bins, normalize=:true, label = "simulation",
title = "Echantillonage pour V(x) = cos($(2 * frequence)*pi*x)")
display(plot!(x, mu, label = "mu"))