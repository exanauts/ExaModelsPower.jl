
using MadNLPHSL, NLPModelsIpopt, NLPModels, LinearAlgebra

cases = [
"pglib_opf_case3_lmbd.m",
"pglib_opf_case5_pjm.m", 
"pglib_opf_case14_ieee.m", 
"pglib_opf_case24_ieee_rts.m", 
"pglib_opf_case30_as.m", 
"pglib_opf_case30_ieee.m", 
"pglib_opf_case5_pjm.m",
"pglib_opf_case14_ieee.m",
"pglib_opf_case24_ieee_rts.m",
"pglib_opf_case30_as.m",
"pglib_opf_case30_ieee.m",
"pglib_opf_case39_epri.m",
"pglib_opf_case57_ieee.m",
"pglib_opf_case60_c.m",
"pglib_opf_case73_ieee_rts.m",
"pglib_opf_case89_pegase.m",
"pglib_opf_case118_ieee.m",
"pglib_opf_case162_ieee_dtc.m",
"pglib_opf_case179_goc.m",
"pglib_opf_case197_snem.m",
"pglib_opf_case200_activ.m",
"pglib_opf_case240_pserc.m",
"pglib_opf_case300_ieee.m",
"pglib_opf_case500_goc.m",
"pglib_opf_case588_sdet.m",
"pglib_opf_case793_goc.m",
"pglib_opf_case1354_pegase.m",
"pglib_opf_case1803_snem.m",
"pglib_opf_case1888_rte.m",
"pglib_opf_case1951_rte.m",
"pglib_opf_case2000_goc.m",
"pglib_opf_case2312_goc.m",
"pglib_opf_case2383wp_k.m",
"pglib_opf_case2736sp_k.m",
"pglib_opf_case2737sop_k.m",
"pglib_opf_case2742_goc.m",
"pglib_opf_case2746wop_k.m",
"pglib_opf_case2746wp_k.m",
"pglib_opf_case2848_rte.m",
"pglib_opf_case2853_sdet.m",
"pglib_opf_case2868_rte.m",
"pglib_opf_case2869_pegase.m",
"pglib_opf_case3012wp_k.m",
"pglib_opf_case3022_goc.m",
"pglib_opf_case3120sp_k.m",
"pglib_opf_case3375wp_k.m",
"pglib_opf_case3970_goc.m",
"pglib_opf_case4020_goc.m",
"pglib_opf_case4601_goc.m",
"pglib_opf_case4619_goc.m",
"pglib_opf_case4661_sdet.m",
"pglib_opf_case4837_goc.m",
"pglib_opf_case4917_goc.m",
"pglib_opf_case5658_epigrids.m",
"pglib_opf_case6468_rte.m",
"pglib_opf_case6470_rte.m",
"pglib_opf_case6495_rte.m",
"pglib_opf_case6515_rte.m",
"pglib_opf_case7336_epigrids.m",
"pglib_opf_case8387_pegase.m",
"pglib_opf_case9241_pegase.m",
"pglib_opf_case9591_goc.m",
"pglib_opf_case10000_goc.m",
"pglib_opf_case10192_epigrids.m",
"pglib_opf_case10480_goc.m",
"pglib_opf_case13659_pegase.m",
"pglib_opf_case19402_goc.m",
"pglib_opf_case20758_epigrids.m",
"pglib_opf_case24464_goc.m",
"pglib_opf_case30000_goc.m",
"pglib_opf_case78484_epigrids.m",]

function termination_code(status::MadNLP.Status)
    if status == MadNLP.SOLVE_SUCCEEDED
        return " "
    elseif status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
        return "a"
    elseif status == MadNLP.DIVERGING_ITERATES || status == MadNLP.DIVERGING_ITERATES
        return "i"
    else
        return "f"
    end
end

function termination_code(status::Symbol)
    if status == :Solve_Succeeded
        return " "
    elseif status == :Solved_To_Acceptable_Level
        return "a"
    elseif status == :Infeasible_Problem_Detected || status == :Diverging_Iterates
        return "i"
    else
        return "f"
    end
end

function evaluate(m, result)
    constraints = similar(result.solution, m.meta.ncon)
    NLPModels.cons!(m, result.solution, constraints)
    return max(
        norm(min.(result.solution .- m.meta.lvar, 0), Inf),
        norm(min.(m.meta.uvar .- result.solution, 0), Inf),
        norm(min.(constraints .- m.meta.lcon, 0), Inf),
        norm(min.(m.meta.ucon .- constraints, 0), Inf)
    )
end

function ipopt_stats(fname)
    output = read(fname, String)
    iter = parse(Int, split(split(output, "Number of Iterations....:")[2], "\n")[1])
    i = parse(Float64,split(split(output, "Total seconds in IPOPT (w/o function evaluations)    =")[2], "\n")[1])
    ad = parse(Float64,split(split(output, "Total seconds in NLP function evaluations            =")[2], "\n")[1])
    tot = i + ad
    return iter, tot, ad
end

using PrettyTables                # at the top of the file
using PrettyTables: tf_latex_booktabs, LatexTableFormat

# after you assemble `methods` and `subs`
function group_boundaries(methods, subs)
    idx  = Int[0, 1]             # columns AFTER which a line is inserted
    col  = 1                 # 1st column = "Case"

    for m in methods
        col += length(subs[m])
        push!(idx, col)      # boundary just after each group
    end

    return idx               # e.g. [4, 7, 9] for groups of size 3,3,2
end




function generate_tex(opf_results; filename="benchmark_results.tex")

    #tf_bt_vlines = deepcopy(tf_latex_booktabs)
    #tf_bt_vlines.col_separator = "|"          # enable vertical bars
    # ————————————————————————————————————————————————————————————
    #  Collect all cases and decide the column layout you want
    # ————————————————————————————————————————————————————————————
    cases = sort(
        [c for c in keys(opf_results) if c != :tol],
        by = x -> parse(Int, match(r"\d+", x).match)
    )

    # explicit ordering ⇣   adjust / reorder as you like
    methods = [:gpu_polar, :gpu_rect, :cpu_polar, :cpu_rect]

    # what sub‑columns belong under each method header?
    subs = Dict(
        :gpu_polar => [:iter, :soltime, :inittime, :adtime,
                       :lintime, :termination, :obj, :cvio],
        :gpu_rect  => [:iter, :soltime, :inittime, :adtime,
                       :lintime, :termination, :obj, :cvio],
        :cpu_polar => [:iter, :soltime,            :adtime,
                       :termination, :obj, :cvio],
        :cpu_rect  => [:iter, :soltime,            :adtime,
                       :termination, :obj, :cvio],
    )

    # ————————————————————————————————————————————————————————————
    #  Build the data matrix (rows) and the two‑level header
    # ————————————————————————————————————————————————————————————
    rows = Any[]
    for case in cases
        row = Any[case]   # left‑hand label
        for m in methods
            for field in subs[m]
                push!(row, get(opf_results[case][m], field, missing))
            end
        end
        push!(rows, row)
    end

    table_data = permutedims(reduce(hcat, rows))   # n_cases × n_cols matrix

    nrows = length(rows)
    hlines = vcat(0, 1, collect(6:5:nrows), nrows+1)

    # two‑level header: top = method name, bottom = metric name
    h_top    = ["Case"]
    h_bottom = [""]

    for m in methods
        n = length(subs[m])                   # how many sub‑columns
        push!(h_top, string(m))               # first column in the span
        append!(h_top, fill("", n-1))    # placeholders for the rest
        append!(h_bottom, string.(subs[m]))   # the sub‑headers themselves
    end

    vlines = group_boundaries(methods, subs)

    # ————————————————————————————————————————————————————————————
    #  Write the .tex file
    # ————————————————————————————————————————————————————————————
    open(filename, "w") do io
        pretty_table(
            io, table_data;
            header = (h_top, h_bottom),
            backend = Val(:latex),
            tf = tf_latex_default,   # ⟹ booktabs rules
            alignment = :c,           # left align everything
            vlines   = vlines,
            hlines   = hlines 
        )
    end
end

function solve_static_cases(cases, tol)

    max_wall_time = Float64(900)

    #Compile time on smallest case
    model_gpu, ~ = opf_model("pglib_opf_case3_lmbd"; backend = CUDABackend())
    ~ = madnlp(model_gpu, tol = tol, max_iter = 3, disable_garbage_collector=true, dual_initialized=true)
    model_gpu, ~ = opf_model("pglib_opf_case3_lmbd"; form = :rect, backend = CUDABackend())
    ~ = madnlp(model_gpu, tol = tol, max_iter = 3, disable_garbage_collector=true, dual_initialized=true)

    model_cpu, ~ = opf_model("pglib_opf_case3_lmbd";)
    ~ = ipopt(model_cpu, tol = tol, max_iter = 3, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma27")
    model_cpu, ~ = opf_model("pglib_opf_case3_lmbd"; form = :rect)
    ~ = ipopt(model_cpu, tol = tol, max_iter = 3, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma27")

    opf_results = Dict()

    opf_results[:tol]=tol

    for case in cases
        case_result = Dict()

        #GPU, Polar
        m_gpu, v_gpu, c_gpu = opf_model(case; backend = CUDABackend())   
        result_gpu = madnlp(m_gpu, tol=tol, max_wall_time = max_wall_time, disable_garbage_collector=true, dual_initialized=true)

        c = evaluate(m_gpu, result_gpu)

        case_result[:gpu_polar] = Dict()
        case_result[:gpu_polar][:iter] = result_gpu.counters.k
        case_result[:gpu_polar][:soltime] = result_gpu.counters.total_time
        case_result[:gpu_polar][:inittime] = result_gpu.counters.init_time
        case_result[:gpu_polar][:adtime] = result_gpu.counters.eval_function_time
        case_result[:gpu_polar][:lintime] = result_gpu.counters.linear_solver_time
        case_result[:gpu_polar][:termination] = termination_code(result_gpu.status)
        case_result[:gpu_polar][:obj] = result_gpu.objective
        case_result[:gpu_polar][:cvio] = c

        #GPU, Rectangular
        m_gpu, v_gpu, c_gpu = opf_model(case; backend = CUDABackend(), form=:rect)   
        result_gpu = madnlp(m_gpu, tol=tol, max_wall_time = max_wall_time, disable_garbage_collector=true, dual_initialized=true)
        
        c = evaluate(m_gpu, result_gpu)

        case_result[:gpu_rect] = Dict()
        case_result[:gpu_rect][:iter] = result_gpu.counters.k
        case_result[:gpu_rect][:soltime] = result_gpu.counters.total_time
        case_result[:gpu_rect][:inittime] = result_gpu.counters.init_time
        case_result[:gpu_rect][:adtime] = result_gpu.counters.eval_function_time
        case_result[:gpu_rect][:lintime] = result_gpu.counters.linear_solver_time
        case_result[:gpu_rect][:termination] = termination_code(result_gpu.status)
        case_result[:gpu_rect][:obj] = result_gpu.objective
        case_result[:gpu_rect][:cvio] = c

        #CPU, Polar
        m_cpu, v_cpu, c_cpu = opf_model(case)
        result_cpu = ipopt(m_cpu, tol = tol, max_wall_time=max_wall_time, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma27", honor_original_bounds = "no", print_timing_statistics = "yes", output_file = "ipopt_output")

        it, tot, ad = ipopt_stats("ipopt_output")

        c = evaluate(m_cpu, result_cpu)
        case_result[:cpu_polar] = Dict()
        case_result[:cpu_polar][:iter] = it
        case_result[:cpu_polar][:soltime] = tot
        case_result[:cpu_polar][:adtime] = ad
        case_result[:cpu_polar][:termination] = termination_code(result_cpu.solver_specific[:internal_msg])
        case_result[:cpu_polar][:obj] = result_cpu.objective
        case_result[:cpu_polar][:cvio] = c

        #CPU, Rectangular
        m_cpu, v_cpu, c_cpu = opf_model(case; form = :rect)
        result_cpu = ipopt(m_cpu, tol = tol, max_wall_time=max_wall_time, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma27", honor_original_bounds = "no", print_timing_statistics = "yes", output_file = "ipopt_output")

        it, tot, ad = ipopt_stats("ipopt_output")

        c = evaluate(m_cpu, result_cpu)
        case_result[:cpu_rect] = Dict()
        case_result[:cpu_rect][:iter] = it
        case_result[:cpu_rect][:soltime] = tot
        case_result[:cpu_rect][:adtime] = ad
        case_result[:cpu_rect][:termination] = termination_code(result_cpu.solver_specific[:internal_msg])
        case_result[:cpu_rect][:obj] = result_cpu.objective
        case_result[:cpu_rect][:cvio] = c

        opf_results[case] = case_result
    end
    
    generate_tex(opf_results; filename = "benchmark_results.tex")

    return opf_results
end