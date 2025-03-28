using Test, ExaModelsPower, MadNLP, MadNLPGPU, KernelAbstractions, CUDA, PowerModels, Ipopt

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

numeric_tolerance(::Type{Float32}) = 1e-2
numeric_tolerance(::Type{Float64}) = 1e-4


function runtests()
    @testset "ExaModelsPower test" begin

        #Test static opf
        #data_pm = PowerModels.parse_file("data/pglib_opf_case3_lmbd.m")

        data, dicts = ExaModelsPower.parse_ac_power_data("pglib_opf_case3_lmbd.m")

        for (T, backend) in CONFIGS
            #Polar tests
            m, v, c = eval(opf_model)("pglib_opf_case3_lmbd.m"; T=T, backend = backend)
            result = madnlp(m; print_level = MadNLP.ERROR)
            va, vm, pg, qg, p, q = v

            nlp_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>result.options.tol, "print_level"=>0)
            result_pm = solve_opf("data/pglib_opf_case3_lmbd.m",ACPPowerModel, nlp_solver)
     
            @testset "$(T), $(backend), polar" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED
                @test isapprox(result.objective, result_pm["objective"], rtol = result.options.tol*100)
                for key in keys(dicts.gen)
                    @test isapprox(Array(solution(result, pg))[dicts.gen[key]], result_pm["solution"]["gen"][string(key)]["pg"], atol = result.options.tol*100)
                    @test isapprox(Array(solution(result, qg))[dicts.gen[key]], result_pm["solution"]["gen"][string(key)]["qg"], atol = result.options.tol*100)
                end
                for key in keys(dicts.bus)
                    @test isapprox(Array(solution(result, va))[dicts.bus[key]], result_pm["solution"]["bus"][string(key)]["va"], atol = result.options.tol*100)
                    @test isapprox(Array(solution(result, vm))[dicts.bus[key]], result_pm["solution"]["bus"][string(key)]["vm"], atol = result.options.tol*100)
                end

                #Branches are encoded differently in solutions, so matches are hard coded
                vars_dict =  Dict("p" => p, "q" => q)
                for st_var in ["p", "q"]
                    var = vars_dict[st_var]
                    @test isapprox(Array(solution(result, var))[1], result_pm["solution"]["branch"]["2"][string(st_var, "f")], atol = result.options.tol*100)
                    @test isapprox(Array(solution(result, var))[2], result_pm["solution"]["branch"]["3"][string(st_var, "f")], atol = result.options.tol*100)
                    @test isapprox(Array(solution(result, var))[3], result_pm["solution"]["branch"]["1"][string(st_var, "f")], atol = result.options.tol*100)
                    @test isapprox(Array(solution(result, var))[4], result_pm["solution"]["branch"]["2"][string(st_var, "t")], atol = result.options.tol*100)
                    @test isapprox(Array(solution(result, var))[5], result_pm["solution"]["branch"]["3"][string(st_var, "t")], atol = result.options.tol*100)
                    @test isapprox(Array(solution(result, var))[6], result_pm["solution"]["branch"]["1"][string(st_var, "t")], atol = result.options.tol*100)
                end
            end

            #Rectangular tests
            m, v, c = eval(opf_model)("pglib_opf_case3_lmbd.m"; T=T, backend = backend, form = :rect)
            result = madnlp(m; print_level = MadNLP.ERROR)
            vr, vim, pg, qg, p, q = v

            nlp_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>result.options.tol, "print_level"=>0)
            result_pm = solve_opf("data/pglib_opf_case3_lmbd.m", ACRPowerModel, nlp_solver)
     
            @testset "$(T), $(backend), rect" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED
                @test isapprox(result.objective, result_pm["objective"], rtol = result.options.tol*100)
                for key in keys(dicts.gen)
                    @test isapprox(Array(solution(result, pg))[dicts.gen[key]], result_pm["solution"]["gen"][string(key)]["pg"], atol = result.options.tol*100)
                    @test isapprox(Array(solution(result, qg))[dicts.gen[key]], result_pm["solution"]["gen"][string(key)]["qg"], atol = result.options.tol*100)
                end
                for key in keys(dicts.bus)
                    @test isapprox(Array(solution(result, vr))[dicts.bus[key]], result_pm["solution"]["bus"][string(key)]["vr"], atol = result.options.tol*100)
                    @test isapprox(Array(solution(result, vim))[dicts.bus[key]], result_pm["solution"]["bus"][string(key)]["vi"], atol = result.options.tol*100)
                end

                #Branches are encoded differently in solutions, so matches are hard coded
                vars_dict =  Dict("p" => p, "q" => q)
                for st_var in ["p", "q"]
                    var = vars_dict[st_var]
                    @test isapprox(Array(solution(result, var))[1], result_pm["solution"]["branch"]["2"][string(st_var, "f")], atol = result.options.tol*100)
                    @test isapprox(Array(solution(result, var))[2], result_pm["solution"]["branch"]["3"][string(st_var, "f")], atol = result.options.tol*100)
                    @test isapprox(Array(solution(result, var))[3], result_pm["solution"]["branch"]["1"][string(st_var, "f")], atol = result.options.tol*100)
                    @test isapprox(Array(solution(result, var))[4], result_pm["solution"]["branch"]["2"][string(st_var, "t")], atol = result.options.tol*100)
                    @test isapprox(Array(solution(result, var))[5], result_pm["solution"]["branch"]["3"][string(st_var, "t")], atol = result.options.tol*100)
                    @test isapprox(Array(solution(result, var))[6], result_pm["solution"]["branch"]["1"][string(st_var, "t")], atol = result.options.tol*100)
                end
            end
        end
    end
end

runtests()
