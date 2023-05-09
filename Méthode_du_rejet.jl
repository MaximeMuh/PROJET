import Plots, Distributions, Random
function rejet(p,q,k,n,b)
    L=[]
    while length(L)!=n
        u=rand()
        x=rand()*b
        if u<=p(x)/(k*q(x))
            push!(L,x)
        end
    end
    return L
end


using Distributions,Plots,QuadGK
b=12
b2=1
p(x)=pdf(Exponential(2),x)
q(x)=pdf(Uniform(0,b),x)
V(x)=exp(-cos(2*pi*x))
q_2(x)=pdf(Uniform(0,b2),x)
k=5
n=10^6
result, error = quadgk(V, 0, 1)
p2(x)=V(x)/result
#recherche de k2
a=[]
for i in 1:10^3
    u=rand()
    push!(a,p2(u)/q_2(u))
end
k2=maximum(a)

L=rejet(p,q,k,n,b)
histogram(L, bins=70,label="Échantillons", normalize=true)
plot!(p, xlims=(0, b), label="Densité cible")

L2=rejet(p2,q_2,k2,n,b2)
histogram(L2, bins=70,label="Échantillons", normalize=true)
plot!(p2, xlims=(0, b2), label="Densité cible")