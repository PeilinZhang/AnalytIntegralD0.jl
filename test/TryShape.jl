#try to calculate from beast mesh data
using BEAST
using CompScienceMeshes
using AnalytIntegralD0
using Test
#an example of mesh: Mesh{3, 3, Float64}(SVector{3, Float64}[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], SVector{3, Int64}[[1, 2, 3]], Dict{SVector{3, Int64}, Int64}())

#make a sphere
Γ = meshsphere(radius=1.0, h=1.0)
#here for testing
vertices = Γ.vertices
faces = Γ.faces
n = length(faces)

#get result from GL
@time begin
    Agl = IntegrateMesh(Γ)
end

#get result from BEAST
@time begin
    X = lagrangecxd0(Γ)
    κ = 0.0
    t = Helmholtz3D.singlelayer(wavenumber=κ) #hypersingular operator(cannot be tested). transpost double layer
    Ab1 = assemble(t,X,X,quadstrat = BEAST.DoubleNumWiltonSauterQStrat(6,7,6,7,10,10,10,10))
end

#checking how close the values are
diff = Agl .- Ab1

# for i in 1:44
#     for j in 1:44
#         println(i, ", ", j)
#         @test diff[i,j] < 1e-7
#     end
# end
for i in 1:44
    for j in 1:44
        if abs(diff[i,j]) > 1e-6
        println("Case index: ", i, ", ", j)
        println("Case vertices x:", vertices[faces[i][1]],vertices[faces[i][2]],vertices[faces[i][3]])
        println("Case vertices y:", vertices[faces[j][1]],vertices[faces[j][2]],vertices[faces[j][3]])
        e1 = Agl[i,j]
        e2 = Ab1[i,j]
        d1 = diff[i,j]
        println("Value from GL: $e1, from BEAST: $e2")
        println("Difference: $d1")
        println()
        end
    end
end
