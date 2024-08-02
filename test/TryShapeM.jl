#test double layer
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
vertices1 = Γ1.vertices
faces1 = Γ1.faces
n1 = length(faces1)
vertices2 = Γ2.vertices
faces2 = Γ2.faces
n2 = length(faces2)

#get result from GL
@time begin
    Agl = IntegrateMesh(Γ1, Γ2, operator = "doublelayer")
end

#get result from BEAST
@time begin
    X1 = lagrangecxd0(Γ1)
    X2 = lagrangecxd0(Γ2)
    κ = 0.0
    t = Helmholtz3D.doublelayer(wavenumber=κ) #hypersingular operator(cannot be tested). transpost double layer
    Ab1 = assemble(t,X1,X2,quadstrat = BEAST.DoubleNumWiltonSauterQStrat(6,7,6,7,10,10,10,10))
end

# #checking how close the values are
# diff1 = Agl .- Ab1
# diff_norm1 = norm(diff1)
# perc_diff1 = norm(diff1)/norm(Ab1)

# @test perc_diff1 < 1e-6
# @test diff_norm1 < 1e-6

# #checking whether result is more accurate
# Ab2 = assemble(t,X1,X2,quadstrat = BEAST.DoubleNumWiltonSauterQStrat(6,7,6,7,30,30,30,30))
# diff2 = Agl .- Ab2
# diff_norm2 = norm(diff2)
# perc_diff2 = norm(diff2)/norm(Ab2)

# # @test perc_diff < 1e-16
# @test diff_norm2 <= diff_norm1

ddd = Agl .- Ab1
for i in 1:n1
    for j in 1:n2
        if abs(ddd[i,j]) > 1e-6 || isnan(Agl[i,j])
        println("Case index: ", i, ", ", j)
        println("Case vertices x:", vertices1[faces1[i][1]],vertices1[faces1[i][2]],vertices1[faces1[i][3]])
        println("Case vertices y:", vertices2[faces2[j][1]],vertices2[faces2[j][2]],vertices2[faces2[j][3]])
        e1 = Agl[i,j]
        e2 = Ab1[i,j]
        d1 = ddd[i,j]
        println("Value from GL: $e1, from BEAST: $e2")
        println("Difference: $d1")
        println()
        end
    end
end