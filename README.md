# ExaModelsPower.jl
ExaModelsPower.jl is an optimal power flow models using ExaModels.jl

## Usage
### Static optimal power flow
```julia
using ExaModelsPower, MadNLP, MadNLPGPU, CUDA

model, vars = opf_model(
    "pglib_opf_case118_ieee.m";
    backend = CUDABackend()
)
result = madnlp(model; tol=1e-6)
```

### Security-constrained optimal power flow
```julia
model, vars = scopf_model(
    "pglib_opf_case118_ieee.m"; contingencies = [1,2],
    backend = CUDABackend()
)
result = madnlp(model; tol=1e-6) # currently failing
```

### Multi-period optimal power flow
```julia
model, vars = mpopf_model(
    "pglib_opf_case118_ieee.m", # static network data
    "/home/sshin/git/ExaModels_Multiperiod/data/case118_onehour_168.Pd", # dynamic load data
    "/home/sshin/git/ExaModels_Multiperiod/data/case118_onehour_168.Qd"; # dynamic load data
    backend = CUDABackend()
)
result = madnlp(model; tol=1e-6)
```
