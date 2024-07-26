# #Calculating using BEAST
using BEAST
using CompScienceMeshes
using LinearAlgebra
using StaticArrays

Γ = meshsphere(radius=1.0, h=1.0)
X = lagrangecxd0(Γ)

κ = 0.0
t = Helmholtz3D.singlelayer(wavenumber=κ)
#hypersingular operator(cannot be tested). transpost double layer

# A = assemble(BEAST.Identity(),X,X)
A = assemble(t,X,X,quadstrat = BEAST.DoubleNumWiltonSauterQStrat(6,7,6,7,10,10,10,10))
# A = assemble(t,X,X)

println(A)