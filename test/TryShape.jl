#try to calculate from beast mesh data
using BEAST
using CompScienceMeshes
using AnalytIntegralD0
#an example of mesh: Mesh{3, 3, Float64}(SVector{3, Float64}[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], SVector{3, Int64}[[1, 2, 3]], Dict{SVector{3, Int64}, Int64}())

#make a sphere
Γ = meshsphere(radius=1.0, h=1.0)
vertices = Γ.vertices
faces = Γ.faces

n = length(faces)
#iterate over faces
A = Matrix{Float64}(undef,n,n)
fill!(A,0)

for i = 1:n
    for j = 1:n
        A[i,j],~,~,~ = GalerkinLaplaceTriGS(vertices[faces[i][1]],vertices[faces[i][2]],vertices[faces[i][3]],vertices[faces[j][1]],vertices[faces[j][2]],vertices[faces[j][3]])
    end
end


