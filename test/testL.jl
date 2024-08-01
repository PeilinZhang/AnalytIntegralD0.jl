#test whether the analytical result La is similar to numerical result Ln; printing out the results for now. Will use isapprox test when everything is debugged.
using AnalytIntegralD0

#define triangles
#try 4, 8
# vertices1 = [ [0.0, -1.0, 0.0],[-0.5844877207446442, -0.562803881114271, -0.584487720744794],[0.0, -0.7071067827963323, -0.7071067795767627] ]
# vertices2 = [ [0.0, -0.7071067827963323, -0.7071067795767627],[0.7071067827963362, -0.2928932185372557, -0.6435942512626279],[0.7071067827963323, -0.7071067795767627, 0.0] ]
# L, M, Ld, Md = GalerkinLaplaceTriGS(vertices1[1], vertices1[2], vertices1[3], vertices2[1], vertices2[2], vertices2[3]) 
# L = L/4/pi; M=M/4/pi
# println("Integral: L: $L, M: $M, Ld: $Ld, Md: $Md")
# println("BEAST_raw: 0.00850919, -0.003990840335676361")

# vertices1 = [ [-0.7071067795767627, -0.7071067827963323, 0.0],[-0.5844877207446442, -0.562803881114271, -0.584487720744794],[0.0, -1.0, 0.0] ]
# vertices2 = [ [0.0, 1.0, 0.0],[0.5844877211148635, 0.5628038803453047, 0.5844877211150131],[0.7071067795767627, 0.7071067827963323, 0.0] ]
# L, M, Ld, Md = GalerkinLaplaceTriGS(vertices1[1], vertices1[2], vertices1[3], vertices2[1], vertices2[2], vertices2[3]) 
# L = L/4/pi; M=M/4/pi
# println("Integral: L: $L, M: $M, Ld: $Ld, Md: $Md")
# println("BEAST: 0.002406664853262945, -0.0012962112782067474")

# vertices1 = [ [-0.6047879252790818, 0.006213853381196598, -0.7963623254918459],[-0.5844877207446442, -0.562803881114271, -0.584487720744794],[-1.0, 0.0, 0.0] ]
# vertices2 = [ [0.0, -0.7071067827963323, -0.7071067795767627],[-0.6047879252790818, 0.006213853381196598, -0.7963623254918459],[0.0, 0.0, -1.0] ]
# L, M, Ld, Md = GalerkinLaplaceTriGS(vertices1[1], vertices1[2], vertices1[3], vertices2[1], vertices2[2], vertices2[3]) 
# L = L/4/pi; M=M/4/pi
# println("Integral: L: $L, M: $M, Ld: $Ld, Md: $Md")
# println("BEAST: x0.002406664853262945x, -0.0032454527906153486")

# vertices1 = [ [0.0, 0.7071067795767627, -0.7071067827963323],[-0.6047879252790818, 0.006213853381196598, -0.7963623254918459],[-0.7071067827963323, 0.7071067795767627, 0.0] ]
# vertices2 = [ [0.0, 0.7071067795767627, -0.7071067827963323],[-0.6047879252790818, 0.006213853381196598, -0.7963623254918459],[-0.7071067827963323, 0.7071067795767627, 0.0] ]
# L, M, Ld, Md = GalerkinLaplaceTriGS(vertices1[1], vertices1[2], vertices1[3], vertices2[1], vertices2[2], vertices2[3]) 
# L = L/4/pi; M=M/4/pi
# println("Integral: L: $L, M: $M, Ld: $Ld, Md: $Md")
# println("BEAST: 0.008674222745422277, -0.005025328957494389")

using LinearAlgebra
using BEAST
using CompScienceMeshes
using LinearAlgebra
using StaticArrays


vertex1 = SVector(0.0, 0.0, 0.0)
vertex2 = SVector(1.0, 0.0, 0.0)
vertex3 = SVector(0.0, 1.0, 0.0)

vertices1 = [vertex1, vertex2, vertex3]
faces1 = [SVector(1, 2, 3)]
mesh1 = Mesh(vertices1, faces1)
X1 = lagrangecxd0(mesh1)

#define a second triangle
vertex12 = SVector(1.0, 1.0, 1.0)
vertex22 = SVector(2.0, 1.0, 1.0)
vertex32 = SVector(1.0, 2.0, 1.0)
vertices2 = [vertex12, vertex22, vertex32]
faces2 = [SVector(1, 2, 3)]
mesh2 = Mesh(vertices2, faces2)
X2 = lagrangecxd0(mesh2)

κ = 0.0
t = Helmholtz3D.singlelayer(wavenumber=κ)

# A = assemble(BEAST.Identity(),X1,X2)
# A = assemble(1,X1,X2)
A = assemble(t,X1,X2,quadstrat = BEAST.DoubleNumWiltonSauterQStrat(6,7,6,7,10,10,10,10))