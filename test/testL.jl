#test whether the analytical result La is similar to numerical result Ln; printing out the results for now. Will use isapprox test when everything is debugged.
using AnalytIntegralD0
using BEAST

# #define triangles
# vertices1 = [ [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0] ]
# vertices2 = [ [1.0, 1.0, 2.0], [2.0, 1.0, 1.0], [1.0, 3.0, 1.0] ]

# n = 10  # number of quadrature points along each triangle edge for numerical calculation
# Ln = integrate_triangles(vertices1, vertices2, n)
# println("Integral: $Ln")


# # #Calculating using BEAST
# using BEAST
# using CompScienceMeshes
# using LinearAlgebra
# using StaticArrays

# vertex1 = SVector(0.0, 0.0, 0.0)
# vertex2 = SVector(1.0, 0.0, 0.0)
# vertex3 = SVector(0.0, 1.0, 0.0)

# vertices = [vertex1, vertex2, vertex3]
# faces = [SVector(1, 2, 3)]
# mesh1 = Mesh(vertices, faces)
# X1 = lagrangecxd0(mesh1)

# #define a second triangle
# vertex1 = SVector(1.0, 1.0, 2.0)
# vertex2 = SVector(2.0, 1.0, 1.0)
# vertex3 = SVector(1.0, 3.0, 1.0)
# vertices = [vertex1, vertex2, vertex3]
# faces = [SVector(1, 2, 3)]
# mesh2 = Mesh(vertices, faces)
# X2 = lagrangecxd0(mesh2)

# κ = 0.0
# t = Helmholtz3D.singlelayer(wavenumber=κ)

# # A = assemble(BEAST.Identity(),X1,X2)
# # A = assemble(1,X1,X2)
# A = assemble(t,X1,X2)

# println(A)
##

#calculate analytical Ln

#define triangles
vertices1 = [ [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0] ]
vertices2 = [ [1.0, 1.0, 1.0], [2.0, 1.0, 1.0], [1.0, 2.0, 0.0] ]
La = GalerkinLaplaceTriGS1(vertices1[1], vertices1[2], vertices1[3], vertices2[1], vertices2[2], vertices2[3]) 
println("Integral: $La")