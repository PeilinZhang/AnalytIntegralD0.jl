function I2t(e,a,h3,h4)
    #Find 2D integral from 2D domian boundaries(1D integrals)
    #There is triangle surface and square surface. This is for the triangle surface
    #the a passed down here should have 2 vectors

    #find parameters h2 and 
    h2, s0 = GSorthogonalization_expan(e,a)

    #Define parmeters of the 1D integrals
    a12 = a[2]
    a14 = a[1]
    a15 = a[1] - a[2]
    e12 = e14 = a[1]*s0[1] + a[2]*s0[2]
    e15 = e12 + a[2]


    #calculate boundary integrals
    I11t = -s0[1] * I1(e12, [a12], h2, h3, h4)
    I12t = -s0[2] * I1(e14, [a14], h2, h3, h4)
    I13t = (1+s0[1]+s0[2])*I1(e15, [a15], h2, h3, h4)

    #calculate 2D integral for square
    I2t = I11t + I12t + I13t

    return I2t 
end