#Zigzag nanoribbon for graphene
using Plots
using LinearAlgebra
chichi = 62 # = 2 Mod 3

G1(kx, t1, t2, a=1) = -2*t2*sin(kx*a*sqrt(3))
G2(kx, t1, t2, a=1) = 2*t1*cos(kx*a*sqrt(3/4))
G3(kx, t1, t2, a=1) = 2*t2*sin(kx*a*sqrt(3/4))
G4(kx, t1, t2, a=1) = t1

function makeH(kx, t1, t2, n=64)

    H1 = zeros(n,n)
    H2 = zeros(n,n)
    H3 = zeros(n,n)
    H4 = zeros(n,n)
    
    for i in 1:n
        
        if i < n-1
            if i%2 == 1
                H1[i,i] = 1
                H2[i,i+1] = 1
                H2[i+1,i] = 1
                H3[i+2,i] = -1
                H3[i,i+2] = 1
            elseif i%2 == 0
                
                H1[i,i] = -1
                H4[i,i+1] = 1
                H4[i+1,i] = 1
                H3[i+2,i] = 1
                H3[i,i+2] = -1
                
            end 
        end 
        
        if i == n-1
            if i%2 == 1
                H1[i,i] = 1
                H2[i,i+1] = 1
                H2[i+1,i] = 1

            elseif i%2 == 0
                
                H1[i,i] = -1
                H4[i,i+1] = 1
                H4[i+1,i] = 1
                
            end   
            
        end
        if i == n 
            if n%2 == 0
                
                H1[i,i] = -1
            end
            if n%2 == 1
                H1[i,i] = 1
            end
        end
    end
    yup = G1(kx,t1,t2) .* H1
    yupp = G2(kx,t1,t2) .* H2
    yuppp = G3(kx,t1,t2) .* H3
    yupppp = G4(kx,t1,t2) .* H4
    #add yup and yuppp to include the NNN terms
    H =  yupp + yupppp# + yup + yuppp
    return H
end

makeH(1,0.5,0.5)

function collecteigen(kx,t1,t2=0.03)
    F = makeH(kx,t1,t2,chichi)
    n = chichi
    yum = []
    CD = eigen(F)
    Vals = CD.values
    for i in 1:n
        push!(yum, round(Vals[i]; sigdigits = 5))
    end
    return real.(yum)
end
function iter_thru(kx,indexx,t1)
    keke = collecteigen(kx,t1)
    return keke[indexx]
end
kx = range(-pi,pi,100)
tt = [0.1,0.3,0.5,0.7]
p1 = plot()
count = 0
for i in 1:chichi
    count += 1
    plot!(kx,iter_thru.(kx,i,1))
end

# display(count)
plot!(title="Zigzag Nanoribbon for graphene",
xlabel="Momentum",ylabel="Energy", legend = false)

