function GalerkinLaplaceTriGS1(x1,x2,x3,y1,y2,y3)

    #1 at the end, meaning it is for single layer potential.
    #Triangle transformation
    # Define the edge vectors
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

    # Calculate the dot product and minimum norm (distance)
    nxy = dot(nx, ny)
    nxy1 = min(norm(nx + ny), norm(nx - ny))

    a1=lx1; a2=-lx3; a3=-ly1; a4=ly3; e4=x1-y1;
    a = [a1, a2, a3, a4]
    h4, s0 = GSorthogonalization_expan(e4, a)
    println()
    println("h4: $h4")
    println()
    e4x = s0[1]*a1 + s0[2]*a2 + s0[3]*a3 + s0[4]*a4

    #Define parmeters of the 3D integrals
    a11 = a12 = a13 = a34 = a[1]
    a21 = a22 = a23 = a36 = a[2]
    a31 = a14 = a15 = a16 = a[3]
    a33 = a24 = a25 = a26 = a[4]
    a32 = a[4] - a[3]
    a35 = a[2] - a[1]
    e31 = e33 = e34 = e36 = a[1]*s0[1] + a[2]*s0[2] + a[3]*s0[3] + a[4]*s0[4]
    e32 = e31 + a[3]
    e35 = e31 + a[1]

    #calculate boundary integrals. Check whether the coefficients equals zero. if zero, then no need calculate lower layers.
    if abs(s0[4]) < zerotol
        I31 = 0
    else
        I31 = -s0[4] * I3(e31, [a11,a21,a31], h4)
    end
    if abs(1 + s0[3] + s0[4]) < zerotol
        I32 = 0
    else
        I32 = (1 + s0[3] + s0[4])*I3(e32, [a12,a22,a32], h4)
    end
    if abs(s0[3]) < zerotol
        I33 = 0
    else
        I33 = -s0[3] * I3(e33, [a14,a24,a34], h4)
    end
    if abs(s0[2]) < zerotol
        I34 = 0
    else
        I34 = -s0[2] * I3(e34, [a14,a24,a34], h4)
    end
    if abs(1 + s0[1] + s0[2]) < zerotol
        I35 = 0
    else
        I35 = (1 + s0[1] + s0[2])*I3(e35, [a15,a25,a35], h4)
    end
    if abs(s0[1]) < zerotol
        I36 = 0
    else
        I36 = -s0[1]*I3(e36, [a16,a26,a36], h4)
    end    

    #calculate 4D integrals
    I4 = I31 + I32 + I33 + I34 + I35 + I36

    L = 4*Ax*Ay*I4

    return L

end