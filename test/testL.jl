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

vertices1 = [ [-0.7071067795767627, -0.7071067827963323, 0.0],[-0.5844877207446442, -0.562803881114271, -0.584487720744794],[0.0, -1.0, 0.0] ]
vertices2 = [ [-0.7071067795767627, -0.7071067827963323, 0.0],[-0.5844877207446442, -0.562803881114271, -0.584487720744794],[0.0, -1.0, 0.0] ]
L, M, Ld, Md = GalerkinLaplaceTriGS(vertices1[1], vertices1[2], vertices1[3], vertices2[1], vertices2[2], vertices2[3]) 
L = L/4/pi; M=M/4/pi
println("Integral: L: $L, M: $M, Ld: $Ld, Md: $Md")
println("BEAST: 0.002406664853262945, -0.0012962112782067474")

# vertices1 = [ [-0.6047879252790818, 0.006213853381196598, -0.7963623254918459],[-0.5844877207446442, -0.562803881114271, -0.584487720744794],[-1.0, 0.0, 0.0] ]
# vertices2 = [ [0.0, -0.7071067827963323, -0.7071067795767627],[-0.6047879252790818, 0.006213853381196598, -0.7963623254918459],[0.0, 0.0, -1.0] ]
# L, M, Ld, Md = GalerkinLaplaceTriGS(vertices1[1], vertices1[2], vertices1[3], vertices2[1], vertices2[2], vertices2[3]) 
# L = L/4/pi; M=M/4/pi
# println("Integral: L: $L, M: $M, Ld: $Ld, Md: $Md")
# println("BEAST: x0.002406664853262945x, -0.0032454527906153486")

# vertices2 = [ [0.0, -1.0, 0.0],[-0.5844877207446442, -0.562803881114271, -0.584487720744794],[0.0, -0.7071067827963323, -0.7071067795767627] ]
# vertices1 = [ [1.566454495553464, 1.046671435041443, -0.6713413590928455],[1.096318417204946, 1.587378564289562, -0.527555894414279],[1.414213559153525, 1.414213565592665, 0.0] ]
# L, M, Ld, Md = GalerkinLaplaceTriGS(vertices1[1], vertices1[2], vertices1[3], vertices2[1], vertices2[2], vertices2[3]) 
# L = L/4/pi; M=M/4/pi
# println("Integral: L: $L, M: $M, Ld: $Ld, Md: $Md")
# println("BEAST: 0.008674222745422277, -0.005025328957494389")