module AnalytIntegralD0

#importing used modules
using LinearAlgebra

#defeine globla constants and variables
const zerotol = 2e-16
case = 0

include("GalerkinLaplaceTriGS1.jl")
include("GalerkinLaplaceTriGS2.jl")
include("GSorthogonalization.jl")
include("cases.jl")
include("I0.jl")
include("I1.jl")
include("I2s.jl")
include("I2t.jl")
include("I3.jl")
include("NumericalTest.jl")

export GSorthogonalization_expan
export GalerkinLaplaceTriGS1 #this is the function that returns single layer integrate_triangles
export GalerkinLaplaceTriGS2
export integrate_triangles #this is for getting a rough numerical result, which will be used for checking whether the analytical result makes sense.

end