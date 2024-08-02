#test single layer
using BEAST
using CompScienceMeshes
using AnalytIntegralD0
using LinearAlgebra
using Test
#an example of mesh: Mesh{3, 3, Float64}(SVector{3, Float64}[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], SVector{3, Int64}[[1, 2, 3]], Dict{SVector{3, Int64}, Int64}())

#make a sphere
Γ1 = meshsphere(radius=1.0, h=1.0)
Γ2 = meshsphere(radius=2.0, h=1.0)
#here for testing

#get result from GL
@time begin
    Agl = IntegrateMesh(Γ1, Γ2)
end

#get result from BEAST
@time begin
    X1 = lagrangecxd0(Γ1)
    X2 = lagrangecxd0(Γ2)
    κ = 0.0
    t = Helmholtz3D.singlelayer(wavenumber=κ) #hypersingular operator(cannot be tested). transpost double layer
    Ab1 = assemble(t,X1,X2,quadstrat = BEAST.DoubleNumWiltonSauterQStrat(6,7,6,7,10,10,10,10))
end

#checking how close the values are
diff1 = Agl .- Ab1
diff_norm1 = norm(diff1)
perc_diff1 = norm(diff1)/norm(Ab1)

@test perc_diff1 < 1e-7
@test diff_norm1 < 1e-7

# checking whether result is more accurate
Ab2 = assemble(t,X1,X2,quadstrat = BEAST.DoubleNumWiltonSauterQStrat(6,7,6,7,30,30,30,30))
diff2 = Agl .- Ab2
diff_norm2 = norm(diff2)
perc_diff2 = norm(diff2)/norm(Ab2)

@test perc_diff2 < 1e-7
@test diff_norm2 <= diff_norm1
