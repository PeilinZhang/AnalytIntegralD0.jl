function GSorthogonalization_expan(e::Vector{Float64},a::Vector{Vector{Float64}})
    # e is a vector; a is a set of vectors that spans a space
    # the returned values:
    # h is the norm of vector component of e that is perpendicular to span(a); 
    # s is the the magnitude of vector components in a such that the component of e that is in span(a) is sum(si0*ai)

    #Orthogonalisation. u is a set of orthogonal vectors such that span(u) = span(a)
    n = length(a)
    u = Vector{Vector{Float64}}(undef, n)
    u[1] = a[1]
    c = zeros(Float64, n, n)
    for i in 2:n
        sum_term = Vector{Float64}(undef, length(u[1]))
        for k in 1:i-1
            ukn = norm(u[k])
            if ukn < zerotol
                ukn = 1
            end
            c[k,i] = dot(u[k],a[i])/(ukn^2)
            sum_term += c[k,i] * u[k]
        end
        u[i] = a[i] - sum_term
    end

    #Find all projection coefficients c
    for i in 1:n
        for j in 1:n
            ujn = norm(u[j])
            if ujn < zerotol
                ujn = 1
            end
            c[j,i] = dot(u[j],a[i])/(ujn^2)
        end
    end

    # return (c,u) #for testing u

    #Find sj0
    b = zeros(Float64, n)
    s = zeros(Float64, n)
    for j in 1:n
        b[j] = dot(u[j],e)
        sum_term = 0.0
        Ajj = c[j,j]*((norm(u[j]))^2)
        if Ajj < zerotol
            s[j] = 0.0 #actually no need this line as this is initialised 0
        else
            for i in j+1:n
                Aji = c[j,i]*((norm(u[j]))^2)
                sum_term += Aji*s[i]
            end
            s[j] = (b[j] - sum_term)/Ajj
        end
    end

    #find h
    sum_term = Vector{Float64}(undef, length(u[1]))
    for i in 1:n
        sum_term += s[i]*a[i]
    end
    h = e - sum_term
    h = norm(h)

    return (h, s)

end