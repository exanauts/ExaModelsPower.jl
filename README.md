# ExaModelsPower.jl
ExaModelsPower.jl is an optimal power flow models using ExaModels.jl

## Usage
```julia
using ExaModelsPower, MadNLP, MadNLPGPU, CUDA

m = opf_model(
	"pglib_opf_case118_ieee.m";
	backend = CUDABackend()
)
madnlp(m; tol=1e-6)

m = scopf_model(
	"pglib_opf_case118_ieee.m";
	backend = CUDABackend()
)
madnlp(m; tol=1e-6)

m = mpopf_model(
	"pglib_opf_case118_ieee.m", # static network data
	"/home/sshin/git/ExaModels_Multiperiod/data/case118_onehour_168.Pd", # dynamic load data
	"/home/sshin/git/ExaModels_Multiperiod/data/case118_onehour_168.Qd"; # dynamic load data
	backend = CUDABackend()
) 
madnlp(m; tol=1e-6)
```
