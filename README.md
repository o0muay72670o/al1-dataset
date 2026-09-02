# AL1 Dataset

This dataset contains the four lowest eigenenergies of a quantum mechanical two-body system, calculated over a grid of AL1 potential parameters and constituent masses using the Gaussian Expansion Method (GEM) implemented in [TwoBody.jl](https://github.com/JuliaFewBody/TwoBody.jl).

## Model

The AL1 potential used here is based on [Meng, Wang, and Oka (2024)](https://doi.org/10.48550/arXiv.2404.01238):

```math
\hat{H}
= \frac{\pmb{p}^2}{2\mu}
+ m_1
+ m_2
- \frac{\kappa}{r}
+ \lambda r
- \Lambda+\frac{2 \pi \kappa^{\prime}}{3 m_1 m_2}
  \frac{\exp(-r^2 / r_0^2)}{\pi^{3 / 2} r_0^2} \boldsymbol{\sigma}_1 \cdot \boldsymbol{\sigma}_2
```

## Solver

The Schrodinger equation was solved using the Gaussian Expansion Method (GEM) implemented in [TwoBody.jl](https://github.com/JuliaFewBody/TwoBody.jl). Twenty Gaussian basis functions were generated in a geometric progression spanning `0.1` to `80.0`. See [`generate.jl`](./generate.jl) for implementation details.

## Dataset Schema

The 52,488 rows in [`data.csv`](./data.csv) form the Cartesian product of the discrete parameter values listed below, with each row containing the four lowest eigenenergies (in GeV) for one parameter combination.

| Column | Unit | Values |
| --- | --- | --- |
| `alpha` | `GeV` | `0.6`, `0.8`, `1.0` |
| `Lambda` | `GeV^2` | `0.10`, `0.15`, `0.20` |
| `sigma` | Dimensionless | `-3`, `1` |
| `Kappa` | Dimensionless | `0.4`, `0.5`, `0.6` |
| `Kappa_` | Dimensionless | `1.80`, `1.85`, `1.90` |
| `A` | `GeV^(B-1)` | `1.60`, `1.65`, `1.70` |
| `B` | Dimensionless | `0.20`, `0.22`, `0.24` |
| `m1`, `m2` | `GeV` | `1.0`, `2.0`, `3.0`, `4.0`, `5.0`, `6.0` |
| `E0` - `E3` | `GeV` | Computed values in `data.csv` |

## Usage

Reading the dataset with Julia:

```julia
using CSV, DataFrames

df = CSV.read("data.csv", DataFrame)
first(df, 5)
```

Reading the dataset with Python:

```python
import pandas as pd

df = pd.read_csv("data.csv")
print(df.head())
```

## Citation

If you use this dataset in your research, please cite it using the metadata provided in [`CITATION.cff`](./CITATION.cff).

## Regenerating

Install [Julia](https://julialang.org/) and the required packages, then run the generation script:

```shell
git clone https://github.com/o0muay72670o/al1-dataset.git
cd al1-dataset
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. generate.jl
```

## Acknowledgments

This work was conducted as part of the [R-IH Summer Program Internship 2026 (Wako Campus)](https://www.riken.jp/en/news_pubs/events/lectures/20260723_2/index.html), hosted by the [Fundamental Physics Data Sharing Unit at the RIKEN Information R&D and Strategy Headquarters (R-IH)](https://r-ih.riken.jp/en/index-2/irdd/physics/).
