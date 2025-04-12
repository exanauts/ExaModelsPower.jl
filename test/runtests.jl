using Test, ExaModelsPower, MadNLP, MadNLPGPU, KernelAbstractions, CUDA, PowerModels, Ipopt, JuMP, ExaModels

include("opf_tests.jl")

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

test_cases = [("data/pglib_opf_case3_lmbd.m", "case3", test_case3),
              ("data/pglib_opf_case5_pjm.m", "case5", test_case5),
              ("data/pglib_opf_case14_ieee.m", "case14", test_case14)]


#Curve = [1, .9, .8, .95, 1]
true_sol_case3_curve = 25384.366465
true_sol_case3_pregen = 29049.351564
true_sol_case5_curve = 78491.04247
true_sol_case5_pregen = 87816.396884

mp_test_cases = [("data/pglib_opf_case3_lmbd.m", "case3", "data/case3_5split.Pd", "data/case3_5split.Qd", true_sol_case3_curve, true_sol_case3_pregen),
                 ("data/pglib_opf_case5_pjm.m", "case5", "data/case5_5split.Pd", "data/case5_5split.Qd", true_sol_case5_curve, true_sol_case5_pregen)]


function example_func(d, srating)
    return d + .2/srating*d^2
end

function runtests()
    @testset "ExaModelsPower test" begin

        for (T, backend) in CONFIGS

            for (filename, case, test_function) in test_cases
                #Test static opf
                data, dicts = ExaModelsPower.parse_ac_power_data(filename)

                #Polar tests
                m, v, c = opf_model(filename; T=T, backend = backend)
                result = madnlp(m; print_level = MadNLP.ERROR)
                va, vm, pg, qg, p, q = v

                nlp_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>result.options.tol, "print_level"=>0)
                result_pm = solve_opf(filename,ACPPowerModel, nlp_solver)
        
                @testset "$(case), static, $(T), $(backend), polar" begin
                    eval(test_function)(result, result_pm, dicts, pg, qg, p, q)
                    test_polar_voltage(result, result_pm, dicts, va, vm)
                end

                #Rectangular tests
                m, v, c = eval(opf_model)(filename; T=T, backend = backend, form = :rect)
                result = madnlp(m; print_level = MadNLP.ERROR)
                vr, vim, pg, qg, p, q = v

                nlp_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>result.options.tol, "print_level"=>0)
                result_pm = solve_opf(filename, ACRPowerModel, nlp_solver)
        
                @testset "$(case), static, $(T), $(backend), rect" begin
                    eval(test_function)(result, result_pm, dicts, pg, qg, p, q)
                    test_rect_voltage(result, result_pm, dicts, vr, vim)
                end
            end
            

            #Test MP
            for (filename, case, Pd_pregen, Qd_pregen, true_sol_curve, true_sol_pregen) in mp_test_cases
                #Curve = [1, .9, .8, .95, 1]

                #Polar
                m, v, c = eval(mpopf_model)(filename, [1, .9, .8, .95, 1]; T = T, backend = backend)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "$(case), MP, $(T), $(backend), curve, polar" begin
                    test_mp_case(result, true_sol_curve)
                end
                #w function
                m, v, c = eval(mpopf_model)(filename, [1, .9, .8, .95, 1], example_func; T = T, backend = backend)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "$(case), MP, $(T), $(backend), curve, polar, func" begin
                    test_mp_case(result, true_sol_curve)
                end
                #Rect
                m, v, c = eval(mpopf_model)(filename, [1, .9, .8, .95, 1]; T = T, backend = backend, form = :rect)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "$(case), MP, $(T), $(backend), curve, rect" begin
                    test_mp_case(result, true_sol_curve)
                end
                #w function
                m, v, c = eval(mpopf_model)(filename, [1, .9, .8, .95, 1], example_func; T = T, backend = backend, form = :rect)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "$(case), MP, $(T), $(backend), curve, rect" begin
                    test_mp_case(result, true_sol_curve)
                end
            

                #Pregenerated Pd and Qd
                #Polar
                m, v, c = eval(mpopf_model)(filename, Pd_pregen, Qd_pregen; T = T, backend = backend)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "$(case), MP, $(T), $(backend), pregen, polar" begin
                    test_mp_case(result, true_sol_pregen)
                end
                #w function
                m, v, c = eval(mpopf_model)(filename, Pd_pregen, Qd_pregen, example_func; T = T, backend = backend)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "$(case), MP, $(T), $(backend), pregen, polar, func" begin
                    test_mp_case(result, true_sol_pregen)
                end

                #Rect
                m, v, c = eval(mpopf_model)(filename, Pd_pregen, Qd_pregen; T = T, backend = backend, form = :rect)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "$(case), MP, $(T), $(backend), pregen, rect" begin
                    test_mp_case(result, true_sol_pregen)
                end
                #w function
                m, v, c = eval(mpopf_model)(filename, Pd_pregen, Qd_pregen, example_func; T = T, backend = backend, form = :rect)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "$(case), MP, $(T), $(backend), pregen, rect, func" begin
                    test_mp_case(result, true_sol_pregen)
                end
            end
            
            #Test MP w storage

            #Curve = [1, .9, .8, .95, 1]
            true_sol = 25358.827525

            #Polar
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd_mod.m", [1, .9, .8, .95, 1]; T = T, backend = backend)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP w storage, $(T), $(backend), curve, polar" begin
                test_mp_case(result, true_sol)
            end
            #Rect
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd_mod.m", [1, .9, .8, .95, 1]; T = T, backend = backend, form = :rect)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP w storage, $(T), $(backend), curve, rect" begin
                test_mp_case(result, true_sol)
            end

            #With function
            true_sol = 25354.331998

            #Polar
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd_mod.m", [1, .9, .8, .95, 1], example_func; T = T, backend = backend)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP w storage, $(T), $(backend), curve, polar, func" begin
                test_mp_case(result, true_sol)
            end
            #Rect
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd_mod.m", [1, .9, .8, .95, 1], example_func; T = T, backend = backend, form = :rect)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP w storage, $(T), $(backend), curve, rect, func" begin
                test_mp_case(result, true_sol)
            end

            #Pregenerated Pd and Qd
            true_sol = 29023.69118

            #Polar
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd_mod.m", "data/case3_5split.Pd", "data/case3_5split.Qd"; T = T, backend = backend)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP w storage, $(T), $(backend), pregen, polar" begin
                test_mp_case(result, true_sol)
            end
            #Rect
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd_mod.m", "data/case3_5split.Pd", "data/case3_5split.Qd"; T = T, backend = backend, form = :rect)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP w storage, $(T), $(backend), pregen, rect" begin
                test_mp_case(result, true_sol)
            end

            #With function
            true_sol = 29019.172473

            #Polar
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd_mod.m", "data/case3_5split.Pd", "data/case3_5split.Qd", example_func; T = T, backend = backend)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP w storage, $(T), $(backend), pregen, polar, func" begin
                test_mp_case(result, true_sol)
            end
            #Rect
            m, v, c = eval(mpopf_model)("pglib_opf_case3_lmbd_mod.m", "data/case3_5split.Pd", "data/case3_5split.Qd", example_func; T = T, backend = backend, form = :rect)
            result = madnlp(m; print_level = MadNLP.ERROR)
            @testset "MP w storage, $(T), $(backend), pregen, rect, func" begin
                test_mp_case(result, true_sol)
            end
        end
    end
end

runtests()
