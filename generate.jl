# https://juliafewbody.github.io/TwoBody.jl/dev/GEM/

using CSV
using DataFrames
using TwoBody

offsets = Float64[]
linear_coefficients = Float64[]
spin_eigenvalues = Float64[]
coulomb_coefficients = Float64[]
spin_spin_coefficients = Float64[]
A = Float64[]
B = Float64[]
results_m1 = Float64[]
results_m2 = Float64[]
results_E0 = Float64[]
results_E1 = Float64[]
results_E2 = Float64[]
results_E3 = Float64[]

for m1 in 1.0:1.0:6.0
for m2 in 1.0:1.0:6.0
for sigma1_dot_sigma2 in [-3, 1]
for Lambda in [0.6, 0.8, 1.0]
for lambda in [0.10, 0.15, 0.20]
for kappa  in [0.4, 0.5, 0.6]
for kappa_prime in [1.80, 1.85, 1.90]
for a in [1.60, 1.65, 1.70]
for b in [0.20, 0.22, 0.24]

	μ = 1 / (1/m1 + 1/m2) 
	r₀ = a*(2*m1*m2/(m1+m2))^(-b)

	H = Hamiltonian(
		RestEnergy(c=1, m=m1),
		RestEnergy(c=1, m=m2),
		Kinetic(hbar=1, m=μ ),
		Coulomb(-kappa),
		Linear(lambda),
		Constant(-Lambda),
		Gaussian(2*π*kappa_prime/3/m1/m2 / ((sqrt(π)*r₀)^3) * sigma1_dot_sigma2, 1/r₀^2),
	)

	BS = GeometricBasisSet(GaussianBasis, 0.1, 80.0, 20)

	res = solve(H, BS)

	push!(offsets,Lambda)
	push!(linear_coefficients, lambda)
	push!(spin_eigenvalues, sigma1_dot_sigma2)
	push!(coulomb_coefficients, kappa)
	push!(spin_spin_coefficients,kappa_prime)
	push!(A, a)
	push!(B, b)
	push!(results_m1, m1)
	push!(results_m2, m2)
	push!(results_E0, res.E[1])
	push!(results_E1, res.E[2])
	push!(results_E2, res.E[3])
	push!(results_E3, res.E[4])

end
end
end
end
end
end
end
end
end

df = DataFrame(Lambda = offsets, lambda = linear_coefficients, sigma1_dot_sigma2 = spin_eigenvalues, kappa = coulomb_coefficients, kappa_prime = spin_spin_coefficients ,A = A, B = B, m1 = results_m1, m2 = results_m2, E0 = results_E0, E1 = results_E1, E2 = results_E2, E3 = results_E3)
println(df)
CSV.write("data.csv", df)
