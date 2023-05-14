# Zigzag nanoribbon for Alpha-T3
# credit to Quantum spin Hall phase transition in the α-T3 lattice, Physical Review B 103, 075419 (2021)
using Plots
using LinearAlgebra
# gr()
chichi = 181 # = 1 Mod 3

function makeH(kx, t1, n, alpha = 0.0)
    
    L1 = 0.0 * t1 # -ve is spin down, +ve is spin up 
    L2 = alpha * L1
    G1(kx, t1, a=1) = -2*t1*cos(kx*a*sqrt(3)/2)
    G2(kx, t1, a=1) = -alpha * t1
    G3(kx, t1, a=1) = -t1
    G4(kx, t1, a=1) = -2*t1*alpha*cos(kx*a*sqrt(3)/2)
    
    n = chichi
    H1 = zeros(ComplexF64,n,n)
    H2 = zeros(ComplexF64,n,n)
    H3 = zeros(ComplexF64,n,n)
    H4 = zeros(ComplexF64,n,n)
    H5 = zeros(ComplexF64,n,n)
    H6 = zeros(ComplexF64,n,n)
    
    for i in 1:n
        if i%3 == 1
            H5[i,i] = ComplexF64(im*(L1-L2)/(3*sqrt(3)))
            if i < n
               H1[i,i+1] = 1 
            end
            if i < n - 1
               H2[i,i+2] = 1
            end
            if i < n-2 
               H6[i,i+3] = ComplexF64(im*(L2-L1)/3*sqrt(3)) 
            end
            if i ==1
                H5[1,1] = ComplexF64(im*L1/3/sqrt(3))
            end
        end
        if i%3 == 0
            H5[i,i] = ComplexF64(im * (L2)/ (3*sqrt(3)))
            if i < n
               H4[i,i+1] = 1 
            end
            if i < n-2
               H6[i,i+3] = ComplexF64(-im*L2/(3*sqrt(3))) 
            end
        end
        if i%3 == 2
            H5[i,i] = ComplexF64(- im * L1/(3*sqrt(3)))
            if i < n-2
               H3[i,i+2] = 1 
            end
            if i < n-2
               H6[i,i+3] = ComplexF64(im*L1/(3*sqrt(3)))
            end
        end
    end
    # Set periodic boundary conditions (Ghost Atom)
    H1[1,n] = 0
    H1[n,1] = 0
    H4[1,n] = 0
    H4[n,1] = 0
  
    M1 = G1(kx,t1) .* H1
    M11 = Hermitian(M1)
    M2 = G2(kx, t1) .* H2
    M22 = Hermitian(M2)
    M3 = G3(kx,t1) .* H3
    M33 = Hermitian(M3)
    M4 = G4(kx,t1) .* H4
    M44 = Hermitian(M4)

    H = M11 + M22 + M33 + M44 
    return H[1:end-1,1:end-1]
end

using Plots

function collect_eigen(kx, t1, n)
    F = makeH(kx, t1, n)
    return eigvals(F)
end

function plot_dispersion(t1, n)
    kx = range(0, pi, 200)
    num_bands = n-1
    
    energy = zeros(length(kx), num_bands)
    
    for i in 1:length(kx)
        eigvals = collect_eigen(kx[i], t1, n)
        energy[i, :] = eigvals[1:num_bands]
    end
    
    plot(kx, energy, xlabel="kx", ylabel="Energy", label="",ylims = (-1.25,1.25))
end

plot_dispersion(1.0, chichi)  
