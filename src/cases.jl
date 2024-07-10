function Identify_cases(h1, h2, h3, h4)
    if h1 > zerotol && h2 < zerotol && h3 < zerotol && h4 < zerotol
        case = 1
    elseif h1 < zerotol && h2 > zerotol && h3 < zerotol && h4 < zerotol
        case = 2
    elseif h1 > zerotol && h2 > zerotol && h3 < zerotol && h4 < zerotol
        case = 3
    elseif h1 > zerotol && h2 < zerotol && h3 > zerotol && h4 < zerotol
        case = 4
    elseif h1 < zerotol && h2 > zerotol && h3 > zerotol && h4 < zerotol
        case = 5
    elseif h1 > zerotol && h2 < zerotol && h3 < zerotol && h4 > zerotol
        case = 6
    elseif h1 < zerotol && h2 > zerotol && h3 < zerotol && h4 > zerotol
        case = 7
    else
        case = 8 #case 8 sould give the error below
    end
    println(case)
    if case >= 8
        error("Error: h does not satisfy any case")
    end
end