#SSH Model
using Plots
using LinearAlgebra
σx=[[0 1]
    [1 0]]
σy=[[0 -im]
    [im 0]]
σz=[[1 0]
    [0 -1]]
# Functions
Rx(k,t1=1,t2=2)=t1-t2*cos(k)
Ry(k,t1=1,t2=2)=-t2*sin(k)

R(k,t1=1,t2=2)=sqrt(Rx(k,t1,t2)^2+Ry(k,t1,t2)^2) #EIGENVALUE#

H(k,t1=1,t2=2)=Rx(k,t1,t2)*σx+Ry(k,t1,t2)*σy #Hamiltonian kernel#
# domain we will calculate on
l=314
ktest = rand()
ks=range(-π,stop=π,length=l)
dk=ks[2]-ks[1]
F = eigen(H(ktest))

print(R((ktest-ktest), 0.67, 0.67))
J=eigen(H((ktest-ktest), 0.67, 0.67))
display(J.vectors)
display(σx)

# Parameters chosen to look at
vas=[i for i in range(0.1, 1, 5)]
# da=[0.0, 0.5,-0.5, 0.5,-0.5, 0.2,-0.2]
va = round.(vas;sigdigits=3)
# w values corresponding to chosen ds
wa=round.(vas/2;sigdigits=3) 

#Energy Dispersion Diagram
plot()
for ii in 1:length(va)
    plot!(ks,R.(ks,va[ii],va[ii]),label="v=$(va[ii]) w=$(va[ii])",legend = true)
    plot!(ks,-R.(ks,va[ii],va[ii]),label="", legend = true)
end
plot!(title="SSH Band diagrams for Varying T1 and constant T2",
xlabel="Momentum",ylabel="Energy", xlims = (-3,3))


#TO PINPOINT THE T that leads to the zero energy state, and to find it's eigenstate.

function findmin(a)
    p = 2
    for i in 1:length(a)
        if a[i] < p
            p = a[i]
        end
    end
    y = findfirst(x -> x == p, a)
    return y
end
        

eiglist = []
i = 0.1
j = 2
k = 10
L = 0.5
va = range(i,j,k)
for ii in 1:length(va)
    eig1 = R(0.0, va[ii], L)
    push!(eiglist,eig1)
end
# display(eiglist)
YES = findmin(eiglist)
Val_T = (YES-1)*((j-1)/10) + i
F = eigen(H(0.0, Val_T, L))
display(F.vectors)


# Parameters chosen to look at
# Energy Dispersion Diagram
kss=[-3,-2,-1,0,1,2,3]
# da=[0.0, 0.5,-0.5, 0.5,-0.5, 0.2,-0.2]
V = [0.1,0.3,0.5,0.8,1,1.5,2]
# w values corresponding to chosen ds
W=V #Floating point error was making it print bad
display(V[5])

plot()
for ii in 1:length(V)
    plot!(kss,R.(kss,V[ii],W[ii]),label="v=$(V[ii]) w=$(W[ii])",legend = true)
    plot!(kss,-R.(kss,V[ii],W[ii]),label="", legend = true)
end
plot!(title="SSH Band diagrams for Varying T1=T2",
xlabel="Momentum",ylabel="Energy")
# display(R.(ks,va[3],wa[3]))
# display(va)
