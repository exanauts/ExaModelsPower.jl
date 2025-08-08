using ExaModelsPower, MadNLP, MadNLPGPU, CUDA, ExaModels
include("benchmark_opf.jl")



#solve_static_cases(sample_cases, 1e-4, "Polar")
#solve_mp_cases(small_sample_cases, curves, 1e-4, "Polar")
#solve_mp_cases(small_sample_cases, curves, 1e-4, "Polar"; storage=true)
solve_stor_cases_comp(small_sample_cases, 1e-4, "Polar")
solve_sc_cases(sc_cases, 1e-2, false)