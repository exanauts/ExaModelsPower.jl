
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
"pglib_opf_case60_c.m",]
#="pglib_opf_case73_ieee_rts.m",
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
"pglib_opf_case78484_epigrids.m",]=#

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

using PrettyTables, Printf                # at the top of the file
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




function generate_tex_opf(opf_results, coords; filename="benchmark_results_opf.tex")
    cases = sort(
        [c for c in keys(opf_results) if c != :tol],
        by = x -> parse(Int, match(r"\d+", x).match)
    )

    methods = ["GPU "*coords, "CPU "*coords]
    subs = Dict(
        "GPU "*coords => [:iter, :soltime, :inittime, :adtime,
                          :lintime, :termination, :obj, :cvio],
        "CPU "*coords => [:iter, :soltime,            :adtime,
                          :termination, :obj, :cvio],
    )

    format_val(field, val) =
    (val === missing || val === nothing) ? missing :
    !(val isa Number) ? string(val) :
    field == :iter ? string(Int(round(val))) :
    field in [:obj, :cvio] ? @sprintf("%.6e", val) :
    @sprintf("%.3e", round(val, sigdigits=4))

    format_k(val) = isnothing(val) || val === missing ? missing : @sprintf("%.1fk", val / 1000)


    rows = Any[]
    for case in cases
        clean_case = replace(case, r"^pglib_opf_case" => "", r"\.m$" => "")

        case_data = opf_results[case]
        nvar = format_k(get(case_data, "nvar", missing))
        ncon = format_k(get(case_data, "ncon", missing))
        row = Any[clean_case, nvar, ncon]
        for m in methods
            for field in subs[m]
                val = get(opf_results[case][m], field, missing)
                push!(row, format_val(field, val))
            end
        end
        push!(rows, row)
    end

    table_data = permutedims(reduce(hcat, rows))  # n_cases × n_cols matrix

    nrows = length(rows)
    hlines = vcat(0, 1, collect(6:5:nrows), nrows+1)

    h_top    = ["Case", "nvars", "ncons"]
    h_bottom = ["",     "",      ""]

    for m in methods
        n = length(subs[m])
        push!(h_top, string(m))
        append!(h_top, fill("", n-1))
        append!(h_bottom, string.(subs[m]))
    end

    function group_boundaries(methods, subs)
        idx = Int[0, 1, 2, 3]
        col = 3
        for m in methods
            col += length(subs[m])
            push!(idx, col)
        end
        return idx
    end

    vlines = group_boundaries(methods, subs)

    open(filename, "w") do io
        pretty_table(
            io, table_data;
            header = (h_top, h_bottom),
            backend = Val(:latex),
            tf = tf_latex_default,
            alignment = :c,
            vlines = vlines,
            hlines = hlines
        )
    end
end

function generate_tex_mpopf(mpopf_results, coords, curve_names; filename="benchmark_results_mpopf.tex")

    # --- Sort cases by number ---
    cases = sort(
        [c for c in keys(mpopf_results) if c != :tol],
        by = x -> parse(Int, match(r"\d+", x).match)
    )

    # --- Build dynamic method names ---
    methods = String[]
    for curve in curve_names
        push!(methods, "GPU $coords $curve")
        push!(methods, "CPU $coords $curve")
    end

    # --- Define what fields each method has ---
    subs = Dict{String, Vector{Symbol}}()
    for method in methods
        if startswith(method, "GPU")
            subs[method] = [:iter, :soltime, :inittime, :adtime, :lintime, :termination, :obj, :cvio]
        elseif startswith(method, "CPU")
            subs[method] = [:iter, :soltime, :adtime, :termination, :obj, :cvio]
        end
    end

    # --- Format values ---
    format_val(field, val) =
        (val === missing || val === nothing) ? missing :
        !(val isa Number) ? string(val) :
        field == :iter ? string(Int(round(val))) :
        field in [:obj, :cvio] ? @sprintf("%.6e", val) :
        @sprintf("%.3e", round(val, sigdigits=4))

    format_k(val) = isnothing(val) || val === missing ? missing : @sprintf("%.1fk", val / 1000)

    # --- Construct rows ---
    rows = Any[]
    for case in cases
        clean_case = replace(case, r"^pglib_opf_case" => "", r"\.m$" => "")

        case_data = mpopf_results[case]
        nvar = format_k(get(case_data, "nvar", missing))
        ncon = format_k(get(case_data, "ncon", missing))
        row = Any[clean_case, nvar, ncon]

        for m in methods
            for field in subs[m]
                val = get(mpopf_results[case][m], field, missing)
                push!(row, format_val(field, val))
            end
        end
        push!(rows, row)
    end

    table_data = permutedims(reduce(hcat, rows))

    # --- Header construction ---
    h_top    = ["Case", "nvars", "ncons"]
    h_bottom = ["",     "",      ""]

    for m in methods
        n = length(subs[m])
        push!(h_top, m)
        append!(h_top, fill("", n-1))
        append!(h_bottom, string.(subs[m]))
    end

    # --- Group boundary vlines ---
    function group_boundaries(methods, subs)
        idx = Int[0, 1, 2, 3]  # vertical lines at case label and after header
        col = 3
        for m in methods
            col += length(subs[m])
            push!(idx, col)
        end
        return idx
    end

    vlines = group_boundaries(methods, subs)

    # --- Horizontal rules every 5 rows ---
    nrows = length(rows)
    hlines = vcat(0, 1, collect(6:5:nrows), nrows+1)

    # --- Output LaTeX ---
    open(filename, "w") do io
        pretty_table(
            io, table_data;
            header = (h_top, h_bottom),
            backend = Val(:latex),
            tf = tf_latex_default,
            alignment = :c,
            vlines = vlines,
            hlines = hlines
        )
    end
end



function solve_static_cases(cases, tol, coords)

    max_wall_time = Float64(900)

    if coords == "Polar"
        form = :polar
    elseif coords == "Rectangular"
        form = :rect
    else
        error("Wrong coords")
    end

    #Compile time on smallest case
    model_gpu, ~ = opf_model("pglib_opf_case3_lmbd"; backend = CUDABackend(), form=form)
    ~ = madnlp(model_gpu, tol = tol, max_iter = 3, disable_garbage_collector=true, dual_initialized=true)

    model_cpu, ~ = opf_model("pglib_opf_case3_lmbd"; form=form)
    ~ = ipopt(model_cpu, tol = tol, max_iter = 3, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma27")
 
    opf_results = Dict()

    opf_results[:tol]=tol

    for case in cases
        case_result = Dict()

        #GPU 
        m_gpu, v_gpu, c_gpu = opf_model(case; backend = CUDABackend(), form=form)   
        result_gpu = madnlp(m_gpu, tol=tol, max_wall_time = max_wall_time, disable_garbage_collector=true, dual_initialized=true)

        c = evaluate(m_gpu, result_gpu)

        case_result["nvar"] = m_gpu.meta.nvar
        case_result["ncon"] = m_gpu.meta.ncon

        case_result["GPU " * coords] = Dict()
        case_result["GPU " * coords][:iter] = result_gpu.counters.k
        case_result["GPU " * coords][:soltime] = result_gpu.counters.total_time
        case_result["GPU " * coords][:inittime] = result_gpu.counters.init_time
        case_result["GPU " * coords][:adtime] = result_gpu.counters.eval_function_time
        case_result["GPU " * coords][:lintime] = result_gpu.counters.linear_solver_time
        case_result["GPU " * coords][:termination] = termination_code(result_gpu.status)
        case_result["GPU " * coords][:obj] = result_gpu.objective
        case_result["GPU " * coords][:cvio] = c

        #CPU
        m_cpu, v_cpu, c_cpu = opf_model(case; form=form)
        result_cpu = ipopt(m_cpu, tol = tol, max_wall_time=max_wall_time, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma27", honor_original_bounds = "no", print_timing_statistics = "yes", output_file = "ipopt_output")

        it, tot, ad = ipopt_stats("ipopt_output")

        c = evaluate(m_cpu, result_cpu)
        case_result["CPU " * coords] = Dict()
        case_result["CPU " * coords][:iter] = it
        case_result["CPU " * coords][:soltime] = tot
        case_result["CPU " * coords][:adtime] = ad
        case_result["CPU " * coords][:termination] = termination_code(result_cpu.solver_specific[:internal_msg])
        case_result["CPU " * coords][:obj] = result_cpu.objective
        case_result["CPU " * coords][:cvio] = c

        opf_results[case] = case_result
    end
    
    generate_tex_opf(opf_results, coords; filename = "benchmark_results_opf.tex")

    return opf_results
end

curves = Dict("Default" => [.64, .60, .58, .56, .58, .62, .66, .73, .81, .88, .95, 1.0, .99, 1.0, 1.0,
    .97, .96, .96, .93, .92, .92, .93, .87, .78, .7],
    "Gentle" => [.88, .90, .88, .86, .87, .88, .9, .92, .93, .95, .97, 1.0, .99, 1.0, 1.0,
    .97, .96, .96, .93, .92, .92, .93, .89, .85, .82],
    "Overburdened" => [.64, .60, .58, .56, .6, .63, .7, .76, .84, .92, .97, 1.01, 1.03, 1.04, 1.06,
    1.08, 1.05, .98, .93, .92, .92, .93, .87, .8, .73])


function solve_mp_cases(cases, curves, tol, coords)

    max_wall_time = Float64(900)

    if coords == "Polar"
        form = :polar
    elseif coords == "Rectangular"
        form = :rect
    else
        error("Wrong coords")
    end

    #Compile time on smallest case
    model_gpu, ~ = mpopf_model("pglib_opf_case3_lmbd", [1,1,1]; backend = CUDABackend(), form=form)
    ~ = madnlp(model_gpu, tol = tol, max_iter = 3, disable_garbage_collector=true, dual_initialized=true)

    model_cpu, ~ = mpopf_model("pglib_opf_case3_lmbd", [1,1,1]; form=form)
    ~ = ipopt(model_cpu, tol = tol, max_iter = 3, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma27")

    mpopf_results = Dict()

    mpopf_results[:tol]=tol

    for case in cases
        println(case)
        case_result = Dict()
        for (curve_name, curve) in curves            

            #GPU
            m_gpu, v_gpu, c_gpu = mpopf_model(case, curve; backend = CUDABackend(), form=form)   
            result_gpu = madnlp(m_gpu, tol=tol, max_wall_time = max_wall_time, disable_garbage_collector=true, dual_initialized=true)

            case_result["nvar"] = m_gpu.meta.nvar
            case_result["ncon"] = m_gpu.meta.ncon

            c = evaluate(m_gpu, result_gpu)
            

            case_result["GPU " * coords * " " *curve_name] = Dict()
            case_result["GPU " * coords * " " *curve_name][:iter] = result_gpu.counters.k
            case_result["GPU " * coords * " " *curve_name][:soltime] = result_gpu.counters.total_time
            case_result["GPU " * coords * " " *curve_name][:inittime] = result_gpu.counters.init_time
            case_result["GPU " * coords * " " *curve_name][:adtime] = result_gpu.counters.eval_function_time
            case_result["GPU " * coords * " " *curve_name][:lintime] = result_gpu.counters.linear_solver_time
            case_result["GPU " * coords * " " *curve_name][:termination] = termination_code(result_gpu.status)
            case_result["GPU " * coords * " " *curve_name][:obj] = result_gpu.objective
            case_result["GPU " * coords * " " *curve_name][:cvio] = c

            #CPU
            m_cpu, v_cpu, c_cpu = mpopf_model(case, curve; form = :rect)
            result_cpu = ipopt(m_cpu, tol = tol, max_wall_time=max_wall_time, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma27", honor_original_bounds = "no", print_timing_statistics = "yes", output_file = "ipopt_output")

            it, tot, ad = ipopt_stats("ipopt_output")

            c = evaluate(m_cpu, result_cpu)
            case_result["CPU " * coords * " " *curve_name] = Dict()
            case_result["CPU " * coords * " " *curve_name][:iter] = it
            case_result["CPU " * coords * " " *curve_name][:soltime] = tot
            case_result["CPU " * coords * " " *curve_name][:adtime] = ad
            case_result["CPU " * coords * " " *curve_name][:termination] = termination_code(result_cpu.solver_specific[:internal_msg])
            case_result["CPU " * coords * " " *curve_name][:obj] = result_cpu.objective
            case_result["CPU " * coords * " " *curve_name][:cvio] = c


            
        end

        mpopf_results[case] = case_result
    end

    curve_names = collect(keys(curves))
    
    generate_tex_mpopf(mpopf_results, coords, curve_names; filename="benchmark_results_mpopf.tex")
    return mpopf_results
end

