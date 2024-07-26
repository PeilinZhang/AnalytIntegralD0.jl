#some unimportant code that I have used. I put them here in case I need them again.

# using BEAST
# using CompScienceMeshes

# Γ = meshsphere(radius=1.0, h=1.0)
# X = lagrangecxd0(Γ)

# κ = 0.0
# t = Helmholtz3D.singlelayer(wavenumber=κ)

# # A = assemble(BEAST.Identity(),X,X)
# A = assemble(t,X,X)


####
# vertex1 = SVector(0.0, 0.0, 0.0)
# vertex2 = SVector(1.0, 0.0, 0.0)
# vertex3 = SVector(0.0, 1.0, 0.0)

# vertices1 = [vertex1, vertex2, vertex3]
# faces1 = [SVector(1, 2, 3)]
# mesh1 = Mesh(vertices1, faces1)
# X1 = lagrangecxd0(mesh1)

# #define a second triangle
# vertex12 = SVector(1.0, 1.0, 1.0)
# vertex22 = SVector(2.0, 1.0, 1.0)
# vertex32 = SVector(1.0, 2.0, 1.0)
# vertices2 = [vertex12, vertex22, vertex32]
# faces2 = [SVector(1, 2, 3)]
# mesh2 = Mesh(vertices2, faces2)
# X2 = lagrangecxd0(mesh2)

# κ = 0.0
# t = Helmholtz3D.singlelayer(wavenumber=κ)

# # A = assemble(BEAST.Identity(),X1,X2)
# # A = assemble(1,X1,X2)
# A = assemble(t,X1,X2)
