using ExaModelsPower, MadNLP, MadNLPGPU, CUDA, ExaModels
println("Using device: ", CUDA.device())
include("benchmark_opf.jl")
solve_mp_cases(cases, curves, 1e-4, "Polar", "GPU")
solve_mp_cases(cases, curves, 1e-4, "Polar", "GPU"; storage = true)
solve_stor_cases_comp(cases, 1e-4, "Polar", "GPU")
solve_sc_cases(sc_cases, 1e-4, false, "GPU")
solve_sc_cases(sc_cases, 1e-4, true, "GPU")