using ExaModelsPower, MadNLP, MadNLPGPU, ExaModels
include("benchmark_opf.jl")
solve_mp_cases(cases, curves, 1e-4, "Polar", "CPU")
solve_mp_cases(cases, curves, 1e-4, "Polar", "CPU"; storage = true)
solve_stor_cases_comp(cases, 1e-4, "Polar", "CPU")
solve_sc_cases(sc_cases, 1e-4, false, "CPU")
solve_sc_cases(sc_cases, 1e-4, true, "CPU")