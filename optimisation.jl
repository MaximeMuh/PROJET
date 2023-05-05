# On importe les bibliotheques
using LinearAlgebra, JuMP, Ipopt, Plots

# On definit les constantes et les fonctions utilisees dans la suite
N = 40 #Nombre de points dans la subdivision de [0,1] (1 exclus)
h = 1 / N # Parametre de la maille dans l'espace d'apporximations des solutions
V(x) = cos(6 * π * x)

# On discretise
X = [i * h for i=0:N-1]
dμ = [exp(-V(x)) for x in X] # On a pas mis la constante de normalisation Z (on en a pas besoin)

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
# On declare que M est symmetrique
M = Symmetric(M)

function construit_A(D)
    """renvoie la matrice A construite a partir du vecteur D"""
    ## Construction de A
    A = zeros(N, N)
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
    # On declare que A est symmetrique
    A = Symmetric(A)

    return A
end

## Fonction objectif
function λ2(D::T...) where {T<:Real}
    """Calcule les valeurs propres generalisees de A et M en fonction du vecteur D"""
    A = construit_A(D)
    return eigen(A, M).values[N-1]
end

function ∇λ2(∇λ::AbstractVector{T}, D::T...) where {T<:Real}
    A = construit_A(D)

    # Vecteur propre associe a λ2
    U = eigen(A, M).vectors[:,N-1]

    # Construction du gradient ∇λ
    for i ∈ 1:N-1
        U_reduit = U[i:i+1]
        ∂diA_reduit = N * dμ[i] * [-1 1;1 -1]
        ∇λ[i] = U_reduit' * ∂diA_reduit * U_reduit
    end
    U_reduit = [U[N], U[1]]
    ∂diA_reduit = N * dμ[N] * [-1 1;1 -1]
    ∇λ[N] = U_reduit' * ∂diA_reduit * U_reduit

    ∇λ /= U' * M * U
    return
end

## Probleme d'optimisation

# Choix du solveur utilise
model = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))

# Declaration des variables
@variable(model, 0.5 <= D[i = 1:N] <= 1)

# Declaration de l'objectif
register(model, :λ2, N, λ2, ∇λ2)
@NLobjective(model, Min, λ2(D...))

optimize!(model)

D_opti = value.(D)
plot(X, D_opti, label = "D optimal")