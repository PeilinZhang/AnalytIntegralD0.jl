module AnalytIntegralD0

#importing used modules
using LinearAlgebra

#defeine globla constants and variables
const zerotol = 1e-100
case = 0

include("GalerkinLaplaceTriGS.jl")
include("GSorthogonalization.jl")
include("cases.jl")
include("I1.jl")
include("I2s.jl")
include("I2t.jl")
include("I3.jl")
include("I4.jl")
include("NumericalTest.jl")

export GSorthogonalization_expan
export GalerkinLaplaceTriGS
export integrate_triangles

end