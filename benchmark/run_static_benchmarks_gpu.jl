using ExaModelsPower, MadNLP, MadNLPGPU, CUDA, ExaModels
println("Using device: ", CUDA.device())
include("benchmark_opf.jl")
for coord in ["Rectangular"]#["Polar", "Rectangular"]
    for tol in [1e-4, 1e-6, 1e-8]
        for style in ["default", "api", "sad"]
            solve_static_cases(cases, tol, coord, "GPU"; case_style = style)
        end
    end
end