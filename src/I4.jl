function I4(e,a)
    #Find 4D integral from 4D domian boundaries(3D integrals)
    #The a passed down here should have 4 vectors

    #find h4 and s0
    h4, s0 = GSorthogonalization_expan(e,a)
    
    #Define parmeters of the 2D integrals
    a11 = a12 = a13 = a34 = a[1]
    a21 = a22 = a23 = a36 = a[2]
    a31 = a14 = a15 = a16 = a[3]
    a32 = a[4] - a[3]
    a35 = a[2] - a[1]
    e31 = e33 = e34 = e36 = a[1]*s0[1] + a[2]*s0[2] + a[3]*s0[3] + a[4]*s0[4]
    e32 = e31 + a[3]
    e35 = e31 + a[1]

    #calculate boundary integrals
    I31 = -s0[4] * I3(e31, [a11,a21,a31], h4)
    I32 = (1 + s0[3] + s0[4])*I3(e32, [a12,a22,a32], h4)
    I33 = -s0[3]*I3(e33, [a14,a24,a34], h4)
    I34 = -s0[2]*I3(e34, [a14,a24,a34], h4)
    I35 = (1 + s0[1] + s0[2])*I3(e35, [a15,a25,a35], h4)
    I36 = -s0[1]*I3(e36, [a16,a26,a36], h4)

    #calculate 3D integrals for 
    I4 = I31 + I32 + I33 + I34 + I35 + I36

    return I4
    
end