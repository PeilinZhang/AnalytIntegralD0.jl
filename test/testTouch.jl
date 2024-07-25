using LinearAlgebra
using AnalytIntegralD0
using Plots

# Generate logarithmically spaced values
ep1min = 10.0
ep1max = 1.0e12
neps = 200
ep1 = range(log10(ep1min), stop=log10(ep1max), length=neps) #There is a problem with using this range function
epsi = 1.0 ./ ep1

# Define points
x1 = [0.0, 0.0, 0.0]
x2 = [1.0, 0.0, 0.0]
x3 = [0.5, sqrt(3)/2, 0.0]
y1 = x2
y2 = x1
y3 = x3

L = zeros(neps)

# Loop over epsilon values
for i in 1:neps
    eps = epsi[i]
    y1[3] = eps
    y2[3] = eps
    L[i], ~, ~, ~ = GalerkinLaplaceTriGS(x1, x2, x3, y1, y2, y3)
    # println(L[i])
end

plot(epsi, L, xscale=:log10, xlabel="1/epsilon", ylabel="L", title="Plot of L against 1/epsilon")
