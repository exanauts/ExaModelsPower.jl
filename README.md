# ExaModelsPower.jl
ExaModelsPower.jl is an optimal power flow models using ExaModels.jl

![CI](https://github.com/exanauts/ExaModelsPower.jl/actions/workflows/ci.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Usage
### Static optimal power flow
```julia
using ExaModelsPower, MadNLP, MadNLPGPU, CUDA

model, vars, cons = opf_model(
    "pglib_opf_case118_ieee.m";
    backend = CUDABackend(),
    form = :polar
)
result = madnlp(model; tol=1e-6)
```

### Security-constrained optimal power flow
```julia
model, vars, cons = scopf_model(
    "pglib_opf_case118_ieee.m"; contingencies = [1,2],
    backend = CUDABackend()
)
result = madnlp(model; tol=1e-6) # currently failing
```

### Multi-period optimal power flow
```julia
model, vars, cons = mpopf_model(
    "pglib_opf_case118_ieee.m", # static network data
    "/home/sshin/git/ExaModels_Multiperiod/data/case118_onehour_168.Pd", # dynamic load data
    "/home/sshin/git/ExaModels_Multiperiod/data/case118_onehour_168.Qd"; # dynamic load data
    backend = CUDABackend()
)
result = madnlp(model; tol=1e-6)

#Alternatively, input a vector to scale baseline demand to generate a demand curve
model, vars, cons = mpopf_model(
    "pglib_opf_case118_ieee.m", # static network data
    [.64, .60, .58, .56, .56, .58, .64, .76, .87, .95, .99, 1.0, .99, 1.0, 1.0,
    .97, .96, .96, .93, .92, .92, .93, .87, .72, .64], #Demand curve
    backend = CUDABackend(),
    corrective_action_ratio = 0.3
)
result = madnlp(model; tol=1e-6)

#mpopf_model can also handle inputs with storage constraints
model, vars, cons = mpopf_model(
    "pglib_opf_case30_ieee_mod.m", # static network data with storage parameters
    "/home/sshin/git/ExaModels_Multiperiod/data/halfhour_30.Pd", # dynamic load data
    "/home/sshin/git/ExaModels_Multiperiod/data/halfhour_30.Qd"; # dynamic load data
    backend = CUDABackend()
)
result = madnlp(model; tol=1e-6)

#Alternatively, provide a smooth function for the charge/discharge efficiency to remove complementarity constraint
function example_func(d, srating)
    return d + .2/srating*d^2
end

model, vars, cons = mpopf_model(
    "pglib_opf_case30_ieee_mod.m", # static network data
    "/home/sshin/git/ExaModels_Multiperiod/data/halfhour_30.Pd", # dynamic load data
    "/home/sshin/git/ExaModels_Multiperiod/data/halfhour_30.Qd"; # dynamic load data
    example_func, #Discharge/charge efficiency modeled along smooth curve
    backend = CUDABackend()
)
result = madnlp(model; tol=1e-6)


#Modified datasets that can be used for testing
#https://github.com/mit-shin-group/multi-period-opf-data
```




