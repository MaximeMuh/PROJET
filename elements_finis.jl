# On importe les bibliotheques
using LinearAlgebra

# On definit les constantes et les fonctions utilisees dans la suite
N = 10 # Nombre de points dans la subdivision de [0,1] (1 exclus)
h = 1 / N # Parametre de la maille dans l'espace d'approximations des solutions
V(x) = cos(2 * π * x)
D_fonction(x) = 1.5 + cos(2 * pi * x)
# On discretise
X = [i * h for i ∈ 0:N-1]
D = D_fonction.(X)
dμ = [exp(-V(x)) for x ∈ X] # On a pas mis la constante de normalisation Z (on en a pas besoin)

## Construction de A
A = zeros(Float64, N, N)
# On remplit la diagonale principale
A[1, 1] = - N * (D[N] * dμ[N] + D[1] * dμ[1])
for i ∈ 1:N-1
    A[i + 1, i + 1] =  - N * (D[i] * dμ[i] + D[i + 1] * dμ[i + 1])
end
# On remplit les diagonales "superieure" et "inferieure"
for i ∈ 0:N-2
    A[i + 1, i + 2] = N * D[i + 1] * dμ[i + 1]
    A[i + 2, i + 1] = N * D[i + 1] * dμ[i + 1]
end
A[1, N] = N * D[N] * dμ[N]
A[N, 1] = N * D[N] * dμ[N]

## Construction de M
M = zeros(Float64, N, N)
# On remplit la diagonale principale
M[1, 1] = 2 * h / 3 * (dμ[N] + dμ[1])
for i ∈ 2:N
    M[i, i] = 2 * h / 3 * (dμ[i - 1] + dμ[i])
end
# On remplit les diagonales "superieure" et "inferieure"
for i ∈ 1:N-1
    M[i, i + 1] = h / 6 * dμ[i]
    M[i + 1, i] = h / 6 * dμ[i]
end
M[1, N] = h / 6 * dμ[N]
M[N, 1] = h / 6 * dμ[N]

# Calcul des valeurs propres generalisees, la premiere valeur propre strictement 
# negative est toujours l'avant derniere valeur renvoyee par eigvals car eigvals
# trie par defaut les valeurs propres par ordre croissant et on sait que 0 est 
# toujours valeur propre d'ordre 1
# On declare que A et M sont symmetriques car cela peut accelerer certains calculs
A = Symmetric(A)
M = Symmetric(M)
λ2 = eigen(A, M).values[N-1]