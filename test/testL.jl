#test whether the analytical result La is similar to numerical result Ln; printing out the results for now. Will use isapprox test when everything is debugged.
using AnalytIntegralD0

#define triangles
vertices1 = [ [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0] ]
vertices2 = [ [1.0, 1.0, 1.0], [2.0, 1.0, 1.0], [1.0, 2.0, 1.0] ]

n = 10  # number of quadrature points along each triangle edge for numerical calculation
Ln = integrate_triangles(vertices1, vertices2, n)
println("Integral: $Ln")

#calculate analytical Ln
La = GalerkinLaplaceTriGS(vertices1[1], vertices1[2], vertices1[3], vertices2[1], vertices2[2], vertices2[3])
println("Integral: $La")