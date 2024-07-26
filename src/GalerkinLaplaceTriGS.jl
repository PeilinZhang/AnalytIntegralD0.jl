function GalerkinLaplaceTriGS(x1,x2,x3,y1,y2,y3)
    #TODO L need to be divided by 4*pi
    #2 at the end, meaning it is for double layer potential.
    #Triangle transformation
    # Define the edge vectors
    x1::Vector{Float64}=x1; x2::Vector{Float64}=x2; x3::Vector{Float64}=x3
    y1::Vector{Float64}=y1; y2::Vector{Float64}=y2; y3::Vector{Float64}=y3
    lx1 = x2 - x1; lx2 = x3 - x2; lx3 = x1 - x3
    ly1 = y2 - y1; ly2 = y3 - y2; ly3 = y1 - y3

    # Calculate the normal vectors and areas
    nx = cross(lx1, -lx3); Ax = norm(nx); nx = nx / Ax; Ax = 0.5 * Ax
    ny = cross(ly1, -ly3); Ay = norm(ny); ny = ny / Ay; Ay = 0.5 * Ay

    # Calculate the lengths of the triangle's edges
    ellx1 = norm(lx1); ellx2 = norm(lx2); ellx3 = norm(lx3)
    elly1 = norm(ly1); elly2 = norm(ly2); elly3 = norm(ly3)

    # Calculate the unit normal vectors
    ncx1 = cross(lx1, nx); ncx2 = cross(lx2, nx); ncx3 = cross(lx3, nx)
    ncx1 = ncx1 / norm(ncx1); ncx2 = ncx2 / norm(ncx2); ncx3 = ncx3 / norm(ncx3)
    ncy1 = cross(ly1, ny); ncy2 = cross(ly2, ny); ncy3 = cross(ly3, ny)
    ncy1 = ncy1 / norm(ncy1); ncy2 = ncy2 / norm(ncy2); ncy3 = ncy3 / norm(ncy3)

    # Calculate the dot product
    nxy = dot(nx, ny)

    a1=lx1; a2=-lx3; a3=-ly1; a4=ly3; e4=x1-y1;
    h4, s0 = GSorthogonalization_expan(e4, [a1, a2, a3, a4])

    #Define parmeters of the 3D integrals
    a11 = a12 = a13 = a34 = a1
    a21 = a22 = a23 = a36 = a2
    a31 = a14 = a15 = a16 = a3
    a33 = a24 = a25 = a26 = a4
    a32 = a4 - a3
    a35 = a2 - a1
    e31 = e33 = e34 = e36 = a1*s0[1] + a2*s0[2] + a3*s0[3] + a4*s0[4]
    e32 = e31 + a3
    e35 = e31 + a1

    #Using h4 to check whether the triangles are parallel
    if h4 < zerotol
        #nondegenerative case
        I31 = I3(e31, [a11,a21,a31], h4)
        I32 = I3(e32, [a12,a22,a32], h4)
        I33 = I3(e33, [a13,a23,a33], h4)
        I34 = I3(e34, [a14,a24,a34], h4)
        I35 = I3(e35, [a15,a25,a35], h4)
        I36 = I3(e36, [a16,a26,a36], h4)

        I4 = - s0[4]*I31 + (1 + s0[3] + s0[4])*I32 - s0[3]*I33 - s0[2]*I34 + (1 + s0[1] + s0[2])*I35 - s0[1]*I36
        L = 4*Ax*Ay*I4 #jacobian is 2A for transformation so here 2Ax * 2Ay = 4 * Ax * Ay

        Fx=6*Ay*(ellx1*ncx1*I34+ellx3*ncx3*I36+ellx2*ncx2*I35) # I'= 3I from the paper formula (4.8), then use (4.5)
        Fy=6*Ax*(elly1*ncy1*I31+elly3*ncy3*I33+elly2*ncy2*I32)

        M = (dot(nx, Fy) - nxy*dot(ny, Fx))/(nxy^2 - 1) #may not work but worked fine till now
        
        # # Below is learned from matlab code. it seems that there is divide by small number problem so a ifelse statement is used for checking whether 1=nxy2 is small
        # nxy2 = nxy * nxy
        # if (1 - nxy2) > 0.5
        #     swi = 0
        #     M = (nxy * dot(ny, Fx) - dot(nx, Fy)) / (1 - nxy2)
        # else
        #     swi = 1
        #     si = sign(nxy)
        #     ep = ny .- si .* nx
        #     eps = norm(ep)
        #     epn = ep ./ eps
        #     Fxy = ((1 - 0.5 * eps * eps) .* Fx .+ Fy) ./ (eps * (1 - 0.25 * eps * eps))
        #     M = si * dot(epn, Fxy)
        # end

    else
        #degenerative case
        delta = -dot(e4, nx)#singed distance between planes
        h4 = abs(delta)

        I32 = I3(e32, [a12,a22,a32], h4)
        I34 = I3(e34, [a14,a24,a34], h4)
        I35 = I3(e35, [a15,a25,a35], h4)
        I36 = I3(e36, [a16,a26,a36], h4)

        I4 = I32 - s0[2]*I34 + (1 + s0[1] + s0[2])*I35 - s0[1]*I36
        L = 4*Ax*Ay*I4

        if h4 < zerotol*norm(e4)
            M = 0
        else
            I32m = I3m(e32, [a12,a22,a32], h4)
            I34m = I3m(e34, [a14,a24,a34], h4)
            I35m = I3m(e35, [a15,a25,a35], h4)
            I36m = I3m(e36, [a16,a26,a36], h4)
    
            I4m = I32m - s0[2]*I34m + (1 + s0[1] + s0[2])*I35m - s0[1]*I36m
            M = 4*Ax*Ay*delta*I4m
        end


        I34d = I3(e4, [a14,a24,a34], 0)
        I35d = I3(e4+a1, [a15,a25,a35], 0)
        I36d = I3(e4, [a16,a26,a36], 0)
        Fx=6*Ay*(ellx1*ncx1*I34d+ellx3*ncx3*I36d+ellx2*ncx2*I35d)
    end

    #finding gradient vector for single layer
    Ld = - Fx - M*nx

    #finding normal derivative for double layer
    H11 = I2s(x1-y1,[lx1, -ly1], 0, 0)
    H12 = I2s(x1-y2,[lx1, -ly2], 0, 0)
    H13 = I2s(x1-y3,[lx1, -ly3], 0, 0)
    H21 = I2s(x2-y1,[lx2, -ly1], 0, 0)
    H22 = I2s(x2-y2,[lx2, -ly2], 0, 0)
    H23 = I2s(x2-y3,[lx2, -ly3], 0, 0)
    H31 = I2s(x3-y1,[lx3, -ly1], 0, 0)
    H32 = I2s(x3-y2,[lx3, -ly2], 0, 0)
    H33 = I2s(x3-y3,[lx3, -ly3], 0, 0)


    Md = -6*(dot(lx1,ly1)*H11 + dot(lx1,ly2)*H12 + dot(lx1,ly3)*H13 + dot(lx2,ly1)*H21 + dot(lx2,ly2)*H22 + dot(lx2,ly3)*H23 + dot(lx3,ly1)*H31 + dot(lx3,ly2)*H32 + dot(lx3,ly3)*H33)

    return L, M, Ld, Md

end