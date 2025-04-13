function test_case3(result, result_pm, result_nlp_pm, dicts, pg, qg, p, q)
    test_static_case(result, result_pm, result_nlp_pm, dicts, pg, qg)

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

function test_static_case(result, result_pm, result_nlp_pm, dicts, pg, qg)
    @test result.status == result_nlp_pm.status
    @test isapprox(result.objective, result_pm["objective"], rtol = result.options.tol*100)
    for key in keys(dicts.gen)
        @test isapprox(Array(solution(result, pg))[dicts.gen[key]], result_pm["solution"]["gen"][string(key)]["pg"], atol = result.options.tol*1000)
        @test isapprox(Array(solution(result, qg))[dicts.gen[key]], result_pm["solution"]["gen"][string(key)]["qg"], atol = result.options.tol*1000)
    end
end

function test_polar_voltage(result, result_pm, dicts, va, vm)
    for key in keys(dicts.bus)
        @test isapprox(Array(solution(result, va))[dicts.bus[key]], result_pm["solution"]["bus"][string(key)]["va"], atol = result.options.tol*100)
        @test isapprox(Array(solution(result, vm))[dicts.bus[key]], result_pm["solution"]["bus"][string(key)]["vm"], rtol = result.options.tol*100)
    end
end

function test_rect_voltage(result, result_pm, dicts, vr, vim)
    for key in keys(dicts.bus)
        @test isapprox(Array(solution(result, vr))[dicts.bus[key]], result_pm["solution"]["bus"][string(key)]["vr"], rtol = result.options.tol*100)
        @test isapprox(Array(solution(result, vim))[dicts.bus[key]], result_pm["solution"]["bus"][string(key)]["vi"], atol = result.options.tol*100)
    end
end

function test_case5(result, result_pm, result_nlp_pm, dicts, pg, qg, p, q)
    test_static_case(result, result_pm, result_nlp_pm, dicts, pg, qg)
end

function test_case14(result, result_pm, result_nlp_pm, dicts, pg, qg, p, q)
    test_static_case(result, result_pm, result_nlp_pm, dicts, pg, qg)
end

function test_float32(m, m64, result, backend)
    x1 = result.solution
    tol = 2.71828^(log(result.options.tol)/2)
    x2 = x1 .* (1 .+ 0.01 .* (2 .* rand(size(x1)) .- 1))
    x3 = x1 .* (1 .+ 0.01 .* (2 .* rand(size(x1)) .- 1))
    for x in [x1, x2, x3]
        @test isapprox(obj(m, x), obj(m64, x), rtol = tol)
        @test isapprox(cons(m, x), cons(m64, x), rtol = tol)
        if backend != CUDABackend()
            @test isapprox(grad(m, x), grad(m64, x), rtol = tol)
            @test isapprox(jac(m, x), jac(m64, x), rtol = tol)
        end
    end
end

function test_mp_case(result, true_sol)
    @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
    @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
end