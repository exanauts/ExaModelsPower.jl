using Test, ExaModelsPower, MadNLP, MadNLPGPU, KernelAbstractions, CUDA, PowerModels, Ipopt, JuMP, ExaModels, NLPModelsJuMP

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

test_cases = [("pglib_opf_case3_lmbd.m", "case3", test_case3),
              ("pglib_opf_case5_pjm.m", "case5", test_case5),
              ("pglib_opf_case14_ieee.m", "case14", test_case14)]

#MP
#Curve = [1, .9, .8, .95, 1]
true_sol_case3_curve = 25384.366465
true_sol_case3_pregen = 29049.351564
true_sol_case5_curve = 78491.04247
true_sol_case5_pregen = 87816.396884
#W storage
true_sol_case3_curve_stor = 25358.827525
true_sol_case3_curve_stor_func = 25354.331998
true_sol_case3_pregen_stor = 29023.69118
true_sol_case3_pregen_stor_func = 29019.172473
true_sol_case5_curve_stor = 68782.0125
true_sol_case5_curve_stor_func = 70235.91846
true_sol_case5_pregen_stor = 79640.08493
true_sol_case5_pregen_stor_func = 79630.14036
mp_test_cases = [("pglib_opf_case3_lmbd.m", "case3", "data/case3_5split.Pd", "data/case3_5split.Qd", true_sol_case3_curve, true_sol_case3_pregen),
                 ("pglib_opf_case5_pjm.m", "case5", "data/case5_5split.Pd", "data/case5_5split.Qd", true_sol_case5_curve, true_sol_case5_pregen)]

mp_stor_test_cases = [("../data/pglib_opf_case3_lmbd_mod.m", "case3", "data/case3_5split.Pd", "data/case3_5split.Qd",
                        true_sol_case3_curve_stor, true_sol_case3_curve_stor_func, true_sol_case3_pregen_stor, true_sol_case3_pregen_stor_func),
                        ("../data/pglib_opf_case5_pjm_mod.m", "case5", "data/case5_5split.Pd", "data/case5_5split.Qd",
                        true_sol_case5_curve_stor, true_sol_case5_curve_stor_func, true_sol_case5_pregen_stor, true_sol_case5_pregen_stor_func)]

function example_func(d, srating)
    return d + .2/srating*d^2
end

function runtests()
    @testset "ExaModelsPower test" begin

        for (T, backend) in CONFIGS

            for (filename, case, test_function) in test_cases
                #Test static opf
                data, dicts = ExaModelsPower.parse_ac_power_data(filename)
                data_pm = PowerModels.parse_file(filename)

                #Polar tests
                m, v, c = opf_model(filename; T=T, backend = backend)
                result = madnlp(m; print_level = MadNLP.ERROR)
                va, vm, pg, qg, p, q = v

                nlp_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>result.options.tol, "print_level"=>0)
                result_pm = solve_opf(filename,ACPPowerModel, nlp_solver)

                
                m_pm = JuMP.Model()
                pm = instantiate_model(data_pm, ACPPowerModel, PowerModels.build_opf, jump_model = m_pm)
                nlp_pm = MathOptNLPModel(m_pm)
                result_nlp_pm = madnlp(nlp_pm; print_level = MadNLP.ERROR)

        
                @testset "$(case), static, $(T), $(backend), polar" begin
                    if T == Float32
                        m64, v64, c64 = opf_model(filename; T=Float64, backend = backend)
                        result = madnlp(m64; print_level = MadNLP.ERROR)
                        test_float32(m, m64, result, backend)
                    else
                        eval(test_function)(result, result_pm, result_nlp_pm, dicts, pg, qg, p, q)
                        test_polar_voltage(result, result_pm, dicts, va, vm)
                    end
                end

                #Rectangular tests
                m, v, c = eval(opf_model)(filename; T=T, backend = backend, form = :rect)
                result = madnlp(m; print_level = MadNLP.ERROR)
                vr, vim, pg, qg, p, q = v

                nlp_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>result.options.tol, "print_level"=>0)
                result_pm = solve_opf(filename, ACRPowerModel, nlp_solver)

                m_pm = JuMP.Model()
                pm = instantiate_model(data_pm, ACRPowerModel, PowerModels.build_opf, jump_model = m_pm)
                nlp_pm = MathOptNLPModel(m_pm)
                result_nlp_pm = madnlp(nlp_pm; print_level = MadNLP.ERROR)
        
                @testset "$(case), static, $(T), $(backend), rect" begin
                    if T == Float32
                        m64, v64, c64 = opf_model(filename; T=Float64, backend = backend, form = :rect)
                        result = madnlp(m64; print_level = MadNLP.ERROR)
                        test_float32(m, m64, result, backend)
                    else
                        eval(test_function)(result, result_pm, result_nlp_pm, dicts, pg, qg, p, q)
                        test_rect_voltage(result, result_pm, dicts, vr, vim)
                    end
                end
            end
            
            #Test MP
            for (filename, case, Pd_pregen, Qd_pregen, true_sol_curve, true_sol_pregen) in mp_test_cases
                #Curve = [1, .9, .8, .95, 1]

                #Polar
                m, v, c = eval(mpopf_model)(filename, [1, .9, .8, .95, 1]; T = T, backend = backend)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "$(case), MP, $(T), $(backend), curve, polar" begin
                    if T == Float32
                        m64, v64, c64 = eval(mpopf_model)(filename, [1, .9, .8, .95, 1]; T=Float64, backend = backend)
                        result = madnlp(m64; print_level = MadNLP.ERROR)
                        test_float32(m, m64, result, backend)
                    else
                        test_mp_case(result, true_sol_curve)
                    end
                end
                #w function
                m, v, c = eval(mpopf_model)(filename, [1, .9, .8, .95, 1], example_func; T = T, backend = backend)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "$(case), MP, $(T), $(backend), curve, polar, func" begin
                    if T == Float32
                        m64, v64, c64 = eval(mpopf_model)(filename, [1, .9, .8, .95, 1], example_func; T=Float64, backend = backend)
                        result = madnlp(m64; print_level = MadNLP.ERROR)
                        test_float32(m, m64, result, backend)
                    else
                        test_mp_case(result, true_sol_curve)
                    end
                end
                #Rect
                m, v, c = eval(mpopf_model)(filename, [1, .9, .8, .95, 1]; T = T, backend = backend, form = :rect)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "$(case), MP, $(T), $(backend), curve, rect" begin
                    if T == Float32
                        m64, v64, c64 = eval(mpopf_model)(filename, [1, .9, .8, .95, 1]; T=Float64, backend = backend, form = :rect)
                        result = madnlp(m64; print_level = MadNLP.ERROR)
                        test_float32(m, m64, result, backend)
                    else
                        test_mp_case(result, true_sol_curve)
                    end
                end
                #w function
                m, v, c = eval(mpopf_model)(filename, [1, .9, .8, .95, 1], example_func; T = T, backend = backend, form = :rect)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "$(case), MP, $(T), $(backend), curve, rect" begin
                    if T == Float32
                        m64, v64, c64 = eval(mpopf_model)(filename, [1, .9, .8, .95, 1], example_func; T=Float64, backend = backend, form = :rect)
                        result = madnlp(m64; print_level = MadNLP.ERROR)
                        test_float32(m, m64, result, backend)
                    else
                        test_mp_case(result, true_sol_curve)
                    end
                end
            

                #Pregenerated Pd and Qd
                #Polar
                m, v, c = eval(mpopf_model)(filename, Pd_pregen, Qd_pregen; T = T, backend = backend)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "$(case), MP, $(T), $(backend), pregen, polar" begin
                    if T == Float32
                        m64, v64, c64 = eval(mpopf_model)(filename, Pd_pregen, Qd_pregen; T=Float64, backend = backend)
                        result = madnlp(m64; print_level = MadNLP.ERROR)
                        test_float32(m, m64, result, backend)
                    else
                        test_mp_case(result, true_sol_pregen)
                    end
                end
                #w function
                m, v, c = eval(mpopf_model)(filename, Pd_pregen, Qd_pregen, example_func; T = T, backend = backend)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "$(case), MP, $(T), $(backend), pregen, polar, func" begin
                    if T == Float32
                        m64, v64, c64 = eval(mpopf_model)(filename, Pd_pregen, Qd_pregen, example_func; T=Float64, backend = backend)
                        result = madnlp(m64; print_level = MadNLP.ERROR)
                        test_float32(m, m64, result, backend)
                    else
                        test_mp_case(result, true_sol_pregen)
                    end
                end

                #Rect
                m, v, c = eval(mpopf_model)(filename, Pd_pregen, Qd_pregen; T = T, backend = backend, form = :rect)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "$(case), MP, $(T), $(backend), pregen, rect" begin
                    if T == Float32
                        m64, v64, c64 = eval(mpopf_model)(filename, Pd_pregen, Qd_pregen; T=Float64, backend = backend, form = :rect)
                        result = madnlp(m64; print_level = MadNLP.ERROR)
                        test_float32(m, m64, result, backend)
                    else
                        test_mp_case(result, true_sol_pregen)
                    end
                end
                #w function
                m, v, c = eval(mpopf_model)(filename, Pd_pregen, Qd_pregen, example_func; T = T, backend = backend, form = :rect)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "$(case), MP, $(T), $(backend), pregen, rect, func" begin
                    if T == Float32
                        m64, v64, c64 = eval(mpopf_model)(filename, Pd_pregen, Qd_pregen, example_func; T=Float64, backend = backend, form = :rect)
                        result = madnlp(m64; print_level = MadNLP.ERROR)
                        test_float32(m, m64, result, backend)
                    else
                        test_mp_case(result, true_sol_pregen)
                    end
                end
            end
            
            #Test MP w storage
            for (filename, case, Pd_pregen, Qd_pregen, true_sol_curve_stor, 
                true_sol_curve_stor_func, true_sol_pregen_stor, true_sol_pregen_stor_func) in mp_stor_test_cases
                
                #Polar
                m, v, c = eval(mpopf_model)(filename, [1, .9, .8, .95, 1]; T = T, backend = backend)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "MP w storage, $(case), $(T), $(backend), curve, polar" begin
                    if T == Float32
                        m64, v64, c64 = eval(mpopf_model)(filename, [1, .9, .8, .95, 1]; T=Float64, backend = backend)
                        result = madnlp(m64; print_level = MadNLP.ERROR)
                        test_float32(m, m64, result, backend)
                    else
                        test_mp_case(result, true_sol_curve_stor)
                    end
                end
                #Rect
                m, v, c = eval(mpopf_model)(filename, [1, .9, .8, .95, 1]; T = T, backend = backend, form = :rect)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "MP w storage, $(case), $(T), $(backend), curve, rect" begin
                    if T == Float32
                        m64, v64, c64 = eval(mpopf_model)(filename, [1, .9, .8, .95, 1]; T=Float64, backend = backend, form = :rect)
                        result = madnlp(m64; print_level = MadNLP.ERROR)
                        test_float32(m, m64, result, backend)
                    else
                        test_mp_case(result, true_sol_curve_stor)
                    end
                end

                #With function
                #Polar
                m, v, c = eval(mpopf_model)(filename, [1, .9, .8, .95, 1], example_func; T = T, backend = backend)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "MP w storage, $(case), $(T), $(backend), curve, polar, func" begin
                    if T == Float32
                        m64, v64, c64 = eval(mpopf_model)(filename, [1, .9, .8, .95, 1], example_func; T=Float64, backend = backend)
                        result = madnlp(m64; print_level = MadNLP.ERROR)
                        test_float32(m, m64, result, backend)
                    else
                        test_mp_case(result, true_sol_curve_stor_func)
                    end
                end
                #Rect
                m, v, c = eval(mpopf_model)(filename, [1, .9, .8, .95, 1], example_func; T = T, backend = backend, form = :rect)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "MP w storage, $(case) $(T), $(backend), curve, rect, func" begin
                    if T == Float32
                        m64, v64, c64 = eval(mpopf_model)(filename, [1, .9, .8, .95, 1], example_func; T=Float64, backend = backend, form = :rect)
                        result = madnlp(m64; print_level = MadNLP.ERROR)
                        test_float32(m, m64, result, backend)
                    else
                        test_mp_case(result, true_sol_curve_stor_func)
                    end
                end

                #Pregenerated Pd and Qd
                #Polar
                m, v, c = eval(mpopf_model)(filename, Pd_pregen, Qd_pregen; T = T, backend = backend)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "MP w storage, $(case), $(T), $(backend), pregen, polar" begin
                    if T == Float32
                        m64, v64, c64 = eval(mpopf_model)(filename, Pd_pregen, Qd_pregen; T=Float64, backend = backend)
                        result = madnlp(m64; print_level = MadNLP.ERROR)
                        test_float32(m, m64, result, backend)
                    else
                        test_mp_case(result, true_sol_pregen_stor)
                    end
                end
                #Rect
                m, v, c = eval(mpopf_model)(filename, Pd_pregen, Qd_pregen; T = T, backend = backend, form = :rect)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "MP w storage, $(case), $(T), $(backend), pregen, rect" begin
                    if T == Float32
                        m64, v64, c64 = eval(mpopf_model)(filename, Pd_pregen, Qd_pregen; T=Float64, backend = backend, form = :rect)
                        result = madnlp(m64; print_level = MadNLP.ERROR)
                        test_float32(m, m64, result, backend)
                    else
                        test_mp_case(result, true_sol_pregen_stor)
                    end
                end

                #With function
                #Polar
                m, v, c = eval(mpopf_model)(filename, Pd_pregen, Qd_pregen, example_func; T = T, backend = backend)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "MP w storage, $(case), $(T), $(backend), pregen, polar, func" begin
                    if T == Float32
                        m64, v64, c64 = eval(mpopf_model)(filename, Pd_pregen, Qd_pregen, example_func; T=Float64, backend = backend)
                        result = madnlp(m64; print_level = MadNLP.ERROR)
                        test_float32(m, m64, result, backend)
                    else
                        test_mp_case(result, true_sol_pregen_stor_func)
                    end
                end
                #Rect
                m, v, c = eval(mpopf_model)(filename, Pd_pregen, Qd_pregen, example_func; T = T, backend = backend, form = :rect)
                result = madnlp(m; print_level = MadNLP.ERROR)
                @testset "MP w storage, $(case), $(T), $(backend), pregen, rect, func" begin
                    if T == Float32
                        m64, v64, c64 = eval(mpopf_model)(filename, Pd_pregen, Qd_pregen, example_func; T=Float64, backend = backend, form = :rect)
                        result = madnlp(m64; print_level = MadNLP.ERROR)
                        test_float32(m, m64, result, backend)
                    else
                        test_mp_case(result, true_sol_pregen_stor_func)
                    end
                end
            end
        end
    end
end

runtests()
