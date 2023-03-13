# On importe les bibliotheques
using LinearAlgebra

# On definit les constantes et les fonctions utilisees dans la suite
N = 1000 # Nombre de points dans la subdivision de [0,1]
h = 1 / N # Parametre de la maille dans l'espace d'apporximations des solutions

D(x) = 1.5 + cos(2 * pi * x)

function dicho(Vp)
    a = 1
    b = N
    nb_iteration_max = UInt64(log2(N)) + 1

    for k=1:nb_iteration_max
        m = (a + b) / 2
        
        if Vp[a] * Vp[m] > 0
            a = m
        else
            b = m
        end

        if Vp[a] * Vp[b] < 0
            return Vp[a]
        end
    end

    return m
end

## Construction de A
A = zeros(Float64, N, N)

# On remplit la diagonale principale
for i = 0:N-1
    A[i + 1, i + 1] = 2 * N * D(i * h)
end

# On remplit les diagonales "superieure" et "inferieure"
for i = 1:N-1
    A[i, i + 1] = - N * D((i-1) * h)
    A[i + 1, i] = - N * D((i-1) * h)
end
A[1, N] = - N * D(0)
A[N, 1] = - N * D(0)


## Construction de M
M = zeros(Float64, N, N)

# On remplit la diagonale principale avec h/2
for i = 1:N
    M[i, i] = h / 2
end

# On remplit les diagonales "superieure" et "inferieure" avec h/4
for i = 1:N-1
    M[i, i + 1] = h / 4
    M[i + 1, i] = h / 4
end
M[1, N] = h / 4
M[N, 1] = h / 4


## Calcul des valeurs propres generalisees
Vp = eigvals(A, M)
 
#Puis on cherche la premiere valeur propre strictement negative par dichotomie
