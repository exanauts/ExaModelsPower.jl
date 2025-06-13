# Introduction

Welcome to the documentation of [ExaModelsPower.jl](https://github.com/exanauts/ExaModelsPower.jl)	
!!! note
    ExaModelsPower runs on julia `VERSION ≥ v"1.9"`
    
!!! warning
	**Please help us improve ExaModelsPower and this documentation!** ExaModelsPower is in the early stage of development, and you may encounter unintended behaviors or missing documentations. If you find anything is not working as intended or documentation is missing, please [open issues](https://github.com/exanauts/ExaModelsPower.jl/issues) or [pull requests](https://github.com/exanauts/ExaModelsPower.jl/pulls) or start [discussions](https://github.com/exanauts/ExaModelsPower.jl/discussions). 
	
## What is ExaModelsPower.jl?
ExaModelsPower.jl is a Julia package for creating optimal power flow (OPF) models. Unlike other OPF modeling frameworks, ExaModelsPower.jl leverages the capabilities of [ExaModels.jl](https://exanauts.github.io/ExaModels.jl/stable/) in order to solve more complex, large-scale versions of the OPF.  ExaModels.jl employs what we call [SIMD](https://en.wikipedia.org/wiki/Single_instruction,_multiple_data) abstraction for [nonlinear programs](https://en.wikipedia.org/wiki/Nonlinear_programming) (NLPs), which allows for the preservation of the parallelizable structure within the model equations, facilitating efficient [automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation) either on the single-thread CPUs, multi-threaded CPUs, as well as [GPU accelerators](https://en.wikipedia.org/wiki/Graphics_processing_unit). More details about SIMD abstraction can be found [here](https://exanauts.github.io/ExaModels.jl/v0.8/simd/). ExaModels.jl compiles (via Julia's compiler) derivative evaluation codes tailored to each computation pattern. Through reverse-mode automatic differentiation using these tailored codes, ExaModels.jl achieves significantly faster derivative evaluation speeds, even when using CPU.

Recent benchmark results demonstrate that derivative evaluation using ExaModels.jl on GPU can be up to two orders of magnitude faster compared to JuMP or AMPL. This enables us to implement more complex versions of the OPF without needing any relaxations. Currently, ExaModelsPower.jl supports developing models for static OPF, multi-period OPF with or without storage, and security constrained OPF. ExaModelsPower.jl also supports a number of flexible options for the user to specify model coordinate system, setup of time-varying demand profiles, and handling of complementarity constraints for storage models. 

## Supported Solvers
ExaModelsPower can be used with any solver that can handle `NLPModel` data type, but several callbacks are not currently implemented, and cause some errors. Currently, it is tested with the following solvers:
- [Ipopt](https://github.com/JuliaSmoothOptimizers/NLPModelsIpopt.jl) (via [NLPModelsIpopt.jl](https://github.com/JuliaSmoothOptimizers/NLPModelsIpopt.jl))
- [MadNLP.jl](https://github.com/MadNLP/MadNLP.jl)

## Documentation Structure
This documentation is structured in the following way.
- This page provides some introductory information about ExaModelsPower.jl
- The step-by-step tutorial of using ExaModelsPower.jl can be found in TK.
- The API Manual provides information on functions provided within ExaModelsPower.jl, as well as information on the constraints and variables implemented in the static and multi-period OPFs


## Supporting ExaModelsPower.jl
- Please report issues and feature requests via the [GitHub issue tracker](https://github.com/exanauts/ExaModelsPower.jl/issues).
- Questions are welcome at [GitHub discussion forum](https://github.com/exanauts/ExaModelsPower.jl/discussions).
