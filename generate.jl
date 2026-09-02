#URL: https://doi.org/10.48550/arXiv.2404.01238
#System:  ηc(1S)
#Gaussian Expansion Method random generator

N_samples = 100
Alpha = Float64[]
Lambda = Float64[]
sigma = Float64[]
Kappa = Float64[]
Kappa_ = Float64[]
A = Float64[]
B = Float64[]
results_m1 = Float64[]
results_m2 = Float64[]
results_E0 = Float64[]
results_E1 = Float64[]
results_E2 = Float64[]
results_E3 = Float64[]

massrange = range(0.1, 6.0 , length=10000)
using CSV
using DataFrames
using TwoBody

for m1 in 1.0:1.0:6.0
for m2 in 1.0:1.0:6.0
for spin in [-3, 1]
for alpha in [0.6, 0.8, 1.0]
for lambda in [0.10, 0.15, 0.20]
for kappa  in [0.4, 0.5, 0.6]
for kappa_ in [1.80, 1.85, 1.90]
for a  in [1.60, 1.65, 1.70]
for b  in [0.20, 0.22, 0.24]

	μ     = 1 / (1/m1 + 1/m2) 
	r₀     = a*(2*m1*m2/(m1+m2))^(-b)

H = Hamiltonian(
	RestEnergy(c=1, m=m1),
	RestEnergy(c=1, m=m2),
	Kinetic(hbar=1, m=μ ),
	Coulomb(-kappa),
	Linear(lambda),
	Constant(-alpha),
	Gaussian(2*π*kappa_/3/m1/m2 / ((sqrt(π)*r₀)^3) * spin, 1/r₀^2),
	)
	BS = GeometricBasisSet(GaussianBasis, 0.1, 80.0, 20)
	res = solve(H, BS)
	push!(Alpha,alpha)
	push!(Lambda, lambda)
	push!(sigma, spin)
	push!(Kappa, kappa)
	push!(Kappa_,kappa_)
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
df = DataFrame(alpha = Alpha, Lambda = Lambda, sigma = sigma, Kappa = Kappa, Kappa_ = Kappa_ ,A = A, B = B, m1 = results_m1, m2 = results_m2, E0 = results_E0, E1 = results_E1, E2 = results_E2, E3 = results_E3)
println(df)
CSV.write("data/$(N_samples)MWOgen.csv", df)
