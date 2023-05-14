#GUANJIHUAN.com EXAMPLE. Since it is a 2x2 Hamiltonian, 4
using Plots
using LinearAlgebra
gr()
σx=[[0 1]
    [1 0]]
σy=[[0 -im]
    [im 0]]
σz=[[1 0]
    [0 -1]]
σI = [[1 0]
    [0 1]]
# Functions
Fz(kx,ky,m=-1,t3=0.5,t2=1) = m + 2*t3*sin(kx)+2*t3*sin(ky)+2*t2*cos(kx+ky)
Fx(kx,ky,t1=1) = 2*t1*cos(kx)
Fy(kx,ky,t1=1) = 2*t1*cos(ky)
H(kx,ky) = Fx(kx,ky).*σx + Fy(kx,ky).*σy + Fz(kx,ky).*σz


#DIFFERENT Mnm function to calculate the summation over m' !=m
# FOR 2x2 MATRIX
vx(kx,ky) = (H(kx+0.00001,ky) - H(kx,ky))/0.00001
vy(kx,ky) = (H(kx,ky+0.00001) - H(kx,ky))/0.00001

function FindEigenVec(kx,ky,i,H)
    F = eigen(H)
    Q = F.vectors
    P = F.values
    
    return Q[:, sortperm(real(P))[i]]
end

function FindEigenVal(kx,ky,i,H)
    F = eigen(H)
    Q = F.values
    New_eigen = sort(real(Q))

    return New_eigen[i]
end

function CalcMnm(kx,ky,index1,index2)
    pun = H(kx,ky)
    xvec = FindEigenVec(kx,ky,index1,pun)
    yvec = FindEigenVec(kx,ky,index2,pun)
    xval = FindEigenVal(kx,ky,index1,pun)
    yval = FindEigenVal(kx,ky,index2,pun)
    total = (transpose(conj(xvec))*(vx(kx,ky))*yvec*transpose(conj(yvec))*(vy(kx,ky))*xvec)/((xval-yval)^2)
    return total
end
    
function getOMEGA(kx,ky,i=1) #Berry 
    a = [1,2]
    deleteat!(a, findall(x->x==i,a))
    p1 = CalcMnm(kx,ky,i,a[1])
    P = -2*im*(p1)
    return P
end

x = pi
n = 200
kx = [i for i in range(-x,x,n)]
ky = [i for i in range(-x,x,n)]

chern = 0

for i in 1:n
    for ii in 1:n
    add = getOMEGA(kx[i],ky[ii])
    chern += add*((2*x/n)^2)
    end
end
chern = chern/(2*pi)

display(chern)
