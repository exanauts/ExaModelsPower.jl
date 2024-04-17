using Test, ExaModelsPower, MadNLP, MadNLPGPU, KernelAbstractions, CUDA

const CONFIGS = [
    (Float32, nothing),
    (Float64, nothing),
    (Float32, CPU()),
    (Float64, CPU()),
]

if CUDA.has_cuda_gpu()
    push!(
        CONFIGS,
        (Float32, CUDABackend()),
    )
    push!(
        CONFIGS,
        (Float64, CUDABackend()),
    )
end

function runtests()
    @testset "ExaModelsExamples test" begin
        for name in ExaModelsPower.NAMES
            for (T, backend) in CONFIGS
                m = eval(name)(; T = T, backend = backend)
                result = madnlp(m; print_level = MadNLP.ERROR)

                @testset "$name" begin
                    @test result.status == MadNLP.SOLVE_SUCCEEDED
                end
            end
        end
    end
end

runtests()
