#test whether the analytical result La is similar to numerical result Ln; printing out the results for now. Will use isapprox test when everything is debugged.
using AnalytIntegralD0
# using BEAST

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
vertices2 = [ [1.0, 1.0, 1.0], [2.0, 1.0, 1.0], [1.0, 2.0, 1.0] ]
L, M, Ld, Md = GalerkinLaplaceTriGS(vertices1[1], vertices1[2], vertices1[3], vertices2[1], vertices2[2], vertices2[3]) 
println("Integral: L: $L, M: $M, Ld: $Ld, Md: $Md")
println("matlab  : 0.1418          0.0474        [-0.0448   -0.0448   -0.0474]         -0.0025")

vertices1 = [ [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0] ]
vertices2 = [ [1.0, 1.0, 1.0], [2.0, 1.0, 1.0], [2.0, 3.0, 2.0] ]
L, M, Ld, Md = GalerkinLaplaceTriGS(vertices1[1], vertices1[2], vertices1[3], vertices2[1], vertices2[2], vertices2[3]) 
println("Integral: L: $L, M: $M, Ld: $Ld, Md: $Md")
println("matlab  : 0.2490          0.0684        [-0.0674   -0.0608   -0.0684]         -0.0189")