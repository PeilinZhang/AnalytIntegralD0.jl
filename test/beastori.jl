using BEAST
using CompScienceMeshes

Γ = meshsphere(radius=1.0, h=1.0)
X = lagrangecxd0(Γ)

κ = 0.0
t = Helmholtz3D.singlelayer(wavenumber=κ)

# A = assemble(BEAST.Identity(),X,X)
A = assemble(t,X,X)
