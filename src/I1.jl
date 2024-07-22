function I1(e1, a1, h2, h3, h4) :: Float64
    #evaluation of 1D integral
    #the a passed down here should have 1 vector in the format of Vector{Vector{Float64}}

    #find parameters h1 and s01
    h1, s01 = GSorthogonalization_expan(e1,a1)
    s01 = s01[1]
    h = R1 = R2 = R4 = 0.0

    #making the smallest two h zero
    h_all = [h1,h2,h3,h4]
    i = sortperm(h_all)
    h_all[i[1]] = 0.0; h_all[i[2]] = 0.0
    h1 = h_all[1]; h2 = h_all[2]; h3 = h_all[3]; h4 = h_all[4]

    #define P
    function P1(s1::Float64) :: Float64
        P1_result = norm(s1+s01)*norm(a1)
        return P1_result
    end

    #define phi
    function phi(index::Int, P_value::Float64; eta = 0.0) :: Float64
        if index == 1
            phi_result = log((P_value+R4)/h)/P_value
        elseif index == 2
            phi_result = atan(eta*P_value, h^2+R4*sqrt(h^2-eta^2))
        elseif index == 3
            phi_result = log((R1+R4)/sqrt(h^2-h1^2))/R1
        elseif index == 4
            phi_result = ((R2*log((R2+R4)/eta)/h2)-log((h2+h)/eta))*eta^2/P_value^2/h2
        end
        if isnan(phi_result) #solving divide by zero
            phi_result = 0.0
        end
        return phi_result
    end

    #Identify cases
    println("h values are: ", h1, " ", h2, " ", h3, " ", h4)
    case = Identify_cases(h1, h2, h3, h4)

    #Define F (3.37)
    function F1(P_value::Float64) :: Float64
        h = sqrt(h1^2+h2^2+h3^2+h4^2)
        R1 = sqrt(P_value^2+h1^2)
        R2 = sqrt(P_value^2+h2^2)
        R4 = sqrt(P_value^2+h^2)
        if case == 1
            F1_result = phi(1,P_value)/6.0
        elseif case == 2
            F1_result = (phi(1,P_value) - 1/(R2+h2))/6.0
        elseif case == 3
            F1_result = (phi(1,P_value) - phi(2,P_value,eta = h1)*h2/h1)/6.0
        elseif case == 4
            F1_result = ((1-h3^2/h1^2)*phi(1, P_value) - 2*h3*phi(2, P_value, eta=h1) + h3^2*phi(3,P_value)/h1^2)/6.0
        elseif case == 5
            F1_result = (h^2*phi(1, P_value)/h2^2 - 1/(R4 + h) - phi(4,P_value,eta=h3))/6.0
        elseif case == 6
            F1_result = ((1-3*h4^2/h1^2)*phi(1,P_value) - (3-h4^2/h1^2)*h4/h1*phi(2,P_value,eta=h1) + 3*h4^2/h1^2*phi(3,P_value) - h4^2/h1^2/(R4+h4))/6.0
        elseif case == 7
            F1_result = ((1-3*h4^2/h2^2)*phi(1, P_value) - (3-h4^2/h1^2)*h4/h1*phi(2,P_value,eta=h1) + 3*h4^2/h1^2*phi(3,P_value) - h4^2/h1^2/(R4+h4))/6.0
        end
        return F1_result
    end
    F11 = F1(P1(1.0))
    F12 = F1(P1(0.0))
    I1 = (1 + s01)*F11 - s01*F12
    println("I1: $I1")

    return I1
end