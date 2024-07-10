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
include("NumericalTest.jl")

export GSorthogonalization_expan
export GalerkinLaplaceTriGS #this is the function that integrates to 4D integral at the top level. The original I4.jl is deleted.
export integrate_triangles #this is for getting a rough numerical result, which will be used for checking whether the analytical result makes sense.

end