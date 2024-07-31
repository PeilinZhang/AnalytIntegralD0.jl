using LinearAlgebra
using Test
using AnalytIntegralD0

zero_tol = 1e-100

# Function to check orthogonality of a set of vectors
function are_orthogonal(vectors)
    n = length(vectors)
    for i in 1:n
        for j in i+1:n
            dot_result = norm(dot(vectors[i], vectors[j]))
            if dot_result > zero_tol
                println(vectors[i])
                println(vectors[j])
                println(dot_result)
                return false
            end
        end
    end
    return true
end

# Function to check if a vector is in the span of a set of vectors
function in_span(v, vectors)
    A = hcat(vectors...)  # Create a matrix with vectors as columns
    rank_A = rank(A)
    rank_Ab = rank(hcat(A, v))
    return rank_A == rank_Ab
end

# Test the 
function test_GSorthogonalization(a)
    e = [0.0,0.0,0.0]
    c, u = GSorthogonalization_expan(e,a)
    println(c)
    println(u)

    # Test for orthogonality
    @test are_orthogonal(u)
    
    # Test that span(u) == span(a)
    for v in a
        @test in_span(v, u)
    end
end

# Run the test
    # Example set of vectors (replace with specific test cases)
a = [
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [1.0, 1.0, 1.0],
    [2.0, 2.0, 0.0]
]

#need to go to the function and discomment the lines that returns u and comments off the lines that returns s and h.
test_GSorthogonalization(a)


# Test GSorthogonalization function
function test_GSorthogonalization(e,a,expected_h,expected_s)
    # expected_h and expected_s are expected outputs

    # Call the function
    h, s = GSorthogonalization_expan(e, a)
    println(h)
    println(s)

    # Check results
    @test isapprox(h, expected_h, atol=1e-6)
    @test all(isapprox.(s, expected_s, atol=1e-6))
end

# Run the test
# # Define test inputs
# e = [1.0, 2.0, 3.0]
# a = [[1.0, 0.0, 0.0],[0.0, 1.0, 0.0], [0.0, 1.0, 0.0],[0.0, 0.0, 0.0]]
# expected_h = 3.0
# expected_s = [1.0, 2.0, 0.0, 0.0]

# test_GSorthogonalization(e,a,expected_h, expected_s)