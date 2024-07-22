function I3(e,a,h4)
    #Find 3D integral from 3D domian boundaries(2D integrals)
    #the 3D shapes are prisms; The a passed down here should have 3 vectors

    #find h3 and s0
    h3, s0 = GSorthogonalization_expan(e,a)
    
    #Define parmeters of the 2D integrals
    a11 = a[1] - a[2]
    a21 = a22 = a23 = a[3]
    a12 = a24 = a25 = a[2]
    a13 = a14 = a15 = a[1]
    e22 = e23 = e25 = a[1]*s0[1]+a[2]*s0[2]+a[3]*s0[3]
    e21 = e22 + a[2]
    e24 = e22 + a[3]

    #calculate boundary integrals
    I21 = (1 + s0[1] + s0[2])*I2s(e21,[a11,a21],h3,h4)
    I22 = -s0[1] * I2s(e22,[a12,a22],h3,h4)
    I23 = -s0[2] * I2s(e23,[a13,a23],h3,h4)
    I24 = (1+s0[3])*I2t(e24,[a14,a24],h3,h4)
    I25 = -s0[3]*I2t(e25,[a15,a25],h3,h4)


    #calculate 3D integrals for 
    I3 = I21 + I22 + I23 + I24 + I25
    println("I3: $I3")

    return I3
    
end