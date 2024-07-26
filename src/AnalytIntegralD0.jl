module AnalytIntegralD0

#importing used modules
using LinearAlgebra

#defeine globla constants and variables
const zerotol = 2e-16
const testMode::Int = 0 #change this into 1 will enable many print parts within the code.
case = 0

include("GalerkinLaplaceTriGS.jl")
include("GSorthogonalization.jl")
include("cases.jl")
include("I0.jl")
include("I1.jl")
include("I2s.jl")
include("I2t.jl")
include("I3.jl")
include("NumericalTest.jl")
include("I0m.jl")
include("I1m.jl")
include("I2sm.jl")
include("I2tm.jl")
include("I3m.jl")

export GSorthogonalization_expan
export GalerkinLaplaceTriGS
export integrate_triangles #this is for getting a rough numerical result, which will be used for checking whether the analytical result makes sense.

end