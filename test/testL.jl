#test whether the analytical result La is similar to numerical result Ln; printing out the results for now. Will use isapprox test when everything is debugged.
using AnalytIntegralD0

#define triangles
# vertices1 = [ [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0] ]
# vertices2 = [ [1.0, 1.0, 1.0], [2.0, 1.0, 1.0], [1.0, 2.0, 1.0] ]
# L, M, Ld, Md = GalerkinLaplaceTriGS(vertices1[1], vertices1[2], vertices1[3], vertices2[1], vertices2[2], vertices2[3]) 
# println("Integral: L: $L, M: $M, Ld: $Ld, Md: $Md")
# println("matlab  : 0.1418          0.0474        [-0.0448   -0.0448   -0.0474]         -0.0025")

vertices1 = [ [-1.0, 0.0, 0.0],[-0.5844877207446442, -0.562803881114271, -0.584487720744794],[-0.7071067795767627, -0.7071067827963323, 0.0] ]
vertices2 = [ [0.7071067795767627, 0.7071067827963323, 0.0],[0.5844877211148635, 0.5628038803453047, 0.5844877211150131],[0.7963623254252495, -0.006213853135373037, 0.6047879253692993] ]
L, M, Ld, Md = GalerkinLaplaceTriGS(vertices1[1], vertices1[2], vertices1[3], vertices2[1], vertices2[2], vertices2[3]) 
L = L/4/pi
println("Integral: L: $L, M: $M, Ld: $Ld, Md: $Md")
println("BEAST4pi: 0.02356776024675801")