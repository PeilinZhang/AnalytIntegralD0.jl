#The file for testing some julia syntax. Just for my learning purpose.
using BEAST
using CompScienceMeshes
using LinearAlgebra
using StaticArrays

vertex1 = SVector(0.0, 0.0, 0.0)
vertex2 = SVector(1.0, 0.0, 0.0)
vertex3 = SVector(0.0, 1.0, 0.0)

vertices = [vertex1, vertex2, vertex3]
faces = [SVector(1, 2, 3)]
mesh1 = Mesh(vertices, faces)
X1 = lagrangecxd0(mesh1)

#define a second triangle
vertex1 = SVector(1.0, 1.0, 2.0)
vertex2 = SVector(2.0, 1.0, 1.0)
vertex3 = SVector(1.0, 3.0, 1.0)
vertices = [vertex1, vertex2, vertex3]
faces = [SVector(1, 2, 3)]
mesh2 = Mesh(vertices, faces)
X2 = lagrangecxd0(mesh2)

κ = 0.0
t = Helmholtz3D.singlelayer(wavenumber=κ)

A = assemble(BEAST.Identity(),X1,X2)
A = assemble(t,X1,X2)