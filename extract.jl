using CSV
using DataFrames
using Random

csv_source = "data.csv"
sample_size = 1000
const seed = 0

function extract_random(csv_source::String, sample_size::Int)
	df = CSV.read(csv_source, DataFrame)
	n_total = nrow(df)
	println("loading $csv_source with $n_total rows")
	rng = MersenneTwister(seed)
	idx = randperm(rng, n_total)[1:sample_size]
	subset = df[idx, :]
    CSV.write("$sample_size data_sample.csv", subset)
end
extract_random(csv_source, sample_size)