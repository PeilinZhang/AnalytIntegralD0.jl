function I2sm(e,a,h3,h4)
    #Find 2D integral from 2D domian boundaries(1D integrals)
    #There is triangle surface and square surface. This is for the sqaure surface

    #find parameters h2 and s0
    h2, s0 = GSorthogonalization_expan(e,a)

    #Define parmeters of the 1D integrals
    a11 = a12 = a[2]
    a13 = a14 = a[1]
    e12 = e14 = a[1]*s0[1] + a[2]*s0[2]
    e11 = e12 + a[1]
    e13 = e12 + a[2]

    #calculate boundary integrals
    if abs(1+s0[1]) < zerotol
        I11s = 0
    else
        I11s = (1+s0[1])*I1m(e11, [a11], h2, h3, h4)
    end
    if abs(s0[1]) < zerotol
        I12s = 0
    else
        I12s = -s0[1] * I1m(e12, [a12], h2, h3, h4)
    end
    if abs(1+s0[2]) < zerotol
        I13s = 0
    else
        I13s = (1+s0[2]) * I1m(e13, [a13], h2, h3, h4)
    end
    if abs(s0[2]) < zerotol
        I14s = 0
    else
        I14s = -s0[2] * I1m(e14, [a14], h2, h3, h4)
    end

    #calculate 2D integral for square
    I2s_sum = I11s + I12s + I13s + I14s
    if testMode == 1
        println("I2m: $I2s_sum")
    end
    return I2s_sum
end