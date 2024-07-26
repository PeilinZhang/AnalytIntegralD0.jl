function IntegrateMesh(Γ::CompScienceMeshes.Mesh{3, 3, Float64})
    #TODO Add operator input: Helmholtz3D etc. Do this when how to integrate into BEAST is figured out because Helmholtz3D is defined in BEAST.
    #Actually I can go into GL and only calculate what is needed to reduce runtime

    #unpack the mesh
    vertices = Γ.vertices
    faces = Γ.faces
    n = length(faces)

    #iterate over faces. TODO output things based on operator
    A = Matrix{Float64}(undef,n,n)
    fill!(A,0)
    for i = 1:n
        for j = 1:n
            # println(vertices[faces[i][1]],vertices[faces[i][2]],vertices[faces[i][3]],vertices[faces[j][1]],vertices[faces[j][2]],vertices[faces[j][3]])
            # if i == 2 && j == 32
            #     error("error is here.")
            # end
            L0,~,~,~ = GalerkinLaplaceTriGS(vertices[faces[i][1]],vertices[faces[i][2]],vertices[faces[i][3]],vertices[faces[j][1]],vertices[faces[j][2]],vertices[faces[j][3]])
            A[i,j] = L0/4/pi
        end
    end

    return A
end