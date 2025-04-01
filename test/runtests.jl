using Test, ExaModelsPower, MadNLP, MadNLPGPU, KernelAbstractions, CUDA, PowerModels, Ipopt, JuMP, ExaModels

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
    @testset "ExaModelsPower test" begin

        for (T, backend) in CONFIGS
            #Test static opf
            data, dicts = ExaModelsPower.parse_ac_power_data("pglib_opf_case3_lmbd.m")

            #Polar tests
            m, v, c = eval(opf_model)("pglib_opf_case3_lmbd.m"; T=T, backend = backend)
            result = madnlp(m; print_level = MadNLP.ERROR)
            va, vm, pg, qg, p, q = v

            nlp_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>result.options.tol, "print_level"=>0)
            result_pm = solve_opf("data/pglib_opf_case3_lmbd.m",ACPPowerModel, nlp_solver)
     
            @testset "static, $(T), $(backend), polar" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
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
     
            @testset "static, $(T), $(backend), rect" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
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
            #Test MP
            #Curve = [1, .9, .8, .95, 1]
            true_sol = 25384.366465

            function example_func(d, srating)
                return d + .2/srating*d^2
            end

            #Polar
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd.m", [1, .9, .8, .95, 1]; T = T, backend = backend)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP, $(T), $(backend), curve, polar" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
                @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
            end
            #w function
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd.m", [1, .9, .8, .95, 1], example_func; T = T, backend = backend)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP, $(T), $(backend), curve, polar, func" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
                @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
            end
            #Rect
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd.m", [1, .9, .8, .95, 1]; T = T, backend = backend, form = :rect)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP, $(T), $(backend), curve, rect" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
                @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
            end
            #w function
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd.m", [1, .9, .8, .95, 1], example_func; T = T, backend = backend, form = :rect)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP, $(T), $(backend), curve, rect" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
                @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
            end

            #Pregenerated Pd and Qd
            true_sol = 29049.351564
            #Polar
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd.m", "data/case3_5split.Pd", "data/case3_5split.Qd"; T = T, backend = backend)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP, $(T), $(backend), pregen, polar" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
                @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
            end
            #w function
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd.m", "data/case3_5split.Pd", "data/case3_5split.Qd", example_func; T = T, backend = backend)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP, $(T), $(backend), pregen, polar, func" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
                @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
            end

            #Rect
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd.m", "data/case3_5split.Pd", "data/case3_5split.Qd"; T = T, backend = backend, form = :rect)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP, $(T), $(backend), pregen, rect" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
                @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
            end
            #w function
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd.m", "data/case3_5split.Pd", "data/case3_5split.Qd", example_func; T = T, backend = backend, form = :rect)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP, $(T), $(backend), pregen, rect, func" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
                @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
            end

            #Test MP w storage

            #Curve = [1, .9, .8, .95, 1]
            true_sol = 25358.827525

            #Polar
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd_mod.m", [1, .9, .8, .95, 1]; T = T, backend = backend)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP w storage, $(T), $(backend), curve, polar" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
                @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
            end
            #Rect
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd_mod.m", [1, .9, .8, .95, 1]; T = T, backend = backend, form = :rect)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP w storage, $(T), $(backend), curve, rect" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
                @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
            end

            #With complementarity constraint
            #Polar
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd_mod.m", [1, .9, .8, .95, 1]; T = T, backend = backend, storage_complementarity_constraint = true)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP w storage, $(T), $(backend), curve, polar, complementarity" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
                @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
            end
            #Rect
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd_mod.m", [1, .9, .8, .95, 1]; T = T, backend = backend, form = :rect, storage_complementarity_constraint = true)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP w storage, $(T), $(backend), curve, rect, complementarity" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
                @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
            end

            #With function
            true_sol = 25354.331998

            #Polar
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd_mod.m", [1, .9, .8, .95, 1], example_func; T = T, backend = backend)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP w storage, $(T), $(backend), curve, polar, func" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
                @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
            end
            #Rect
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd_mod.m", [1, .9, .8, .95, 1], example_func; T = T, backend = backend, form = :rect)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP w storage, $(T), $(backend), curve, rect, func" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
                @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
            end

            #Pregenerated Pd and Qd
            true_sol = 29023.69118

            #Polar
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd_mod.m", "data/case3_5split.Pd", "data/case3_5split.Qd"; T = T, backend = backend)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP w storage, $(T), $(backend), pregen, polar" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
                @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
            end
            #Rect
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd_mod.m", "data/case3_5split.Pd", "data/case3_5split.Qd"; T = T, backend = backend, form = :rect)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP w storage, $(T), $(backend), pregen, rect" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
                @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
            end

            #With complementarity
            #Polar
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd_mod.m", "data/case3_5split.Pd", "data/case3_5split.Qd"; T = T, backend = backend, storage_complementarity_constraint = true)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP w storage, $(T), $(backend), pregen, polar, complementarity" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
                @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
            end
            #Rect
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd_mod.m", "data/case3_5split.Pd", "data/case3_5split.Qd"; T = T, backend = backend, form = :rect, storage_complementarity_constraint = true)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP w storage, $(T), $(backend), pregen, rect, complementarity" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
                @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
            end

            #With function
            true_sol = 29019.172473

            #Polar
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd_mod.m", "data/case3_5split.Pd", "data/case3_5split.Qd", example_func; T = T, backend = backend)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP w storage, $(T), $(backend), pregen, polar, func" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
                @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
            end
            #Rect
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd_mod.m", "data/case3_5split.Pd", "data/case3_5split.Qd", example_func; T = T, backend = backend, form = :rect)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP w storage, $(T), $(backend), pregen, rect, func" begin
                @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
                @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
            end
        end
    end
end

runtests()
