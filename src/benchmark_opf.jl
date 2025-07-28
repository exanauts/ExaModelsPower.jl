
using MadNLPHSL, NLPModelsIpopt, NLPModels, LinearAlgebra, CSV, DataFrames, PrettyTables, Printf, Plots, SolverBenchmark, MadNCL, HybridKKT
using PrettyTables: tf_latex_booktabs, LatexTableFormat


cases = [
"pglib_opf_case3_lmbd",
"pglib_opf_case5_pjm", 
"pglib_opf_case14_ieee", 
"pglib_opf_case24_ieee_rts", 
"pglib_opf_case30_as", 
"pglib_opf_case30_ieee", 
"pglib_opf_case5_pjm",
"pglib_opf_case14_ieee",]
#="pglib_opf_case24_ieee_rts",
"pglib_opf_case30_as",
"pglib_opf_case30_ieee",
"pglib_opf_case39_epri",
"pglib_opf_case57_ieee",
"pglib_opf_case60_c",
"pglib_opf_case73_ieee_rts",
"pglib_opf_case89_pegase",
"pglib_opf_case118_ieee",
"pglib_opf_case162_ieee_dtc",
"pglib_opf_case179_goc",
"pglib_opf_case197_snem",
"pglib_opf_case200_activ",
"pglib_opf_case240_pserc",
"pglib_opf_case300_ieee",
"pglib_opf_case500_goc",
"pglib_opf_case588_sdet",
"pglib_opf_case793_goc",
"pglib_opf_case1354_pegase",
"pglib_opf_case1803_snem",
"pglib_opf_case1888_rte",
"pglib_opf_case1951_rte",
"pglib_opf_case2000_goc",
"pglib_opf_case2312_goc",
"pglib_opf_case2383wp_k",
"pglib_opf_case2736sp_k",
"pglib_opf_case2737sop_k",
"pglib_opf_case2742_goc",
"pglib_opf_case2746wop_k",
"pglib_opf_case2746wp_k",
"pglib_opf_case2848_rte",
"pglib_opf_case2853_sdet",
"pglib_opf_case2868_rte",
"pglib_opf_case2869_pegase",
"pglib_opf_case3012wp_k",
"pglib_opf_case3022_goc",
"pglib_opf_case3120sp_k",
"pglib_opf_case3375wp_k",
"pglib_opf_case3970_goc",
"pglib_opf_case4020_goc",
"pglib_opf_case4601_goc",
"pglib_opf_case4619_goc",
"pglib_opf_case4661_sdet",
"pglib_opf_case4837_goc",
"pglib_opf_case4917_goc",
"pglib_opf_case5658_epigrids",
"pglib_opf_case6468_rte",
"pglib_opf_case6470_rte",
"pglib_opf_case6495_rte",
"pglib_opf_case6515_rte",
"pglib_opf_case7336_epigrids",
"pglib_opf_case8387_pegase",
"pglib_opf_case9241_pegase",
"pglib_opf_case9591_goc",
"pglib_opf_case10000_goc",
"pglib_opf_case10192_epigrids",
"pglib_opf_case10480_goc",
"pglib_opf_case13659_pegase",
"pglib_opf_case19402_goc",
"pglib_opf_case20758_epigrids",
"pglib_opf_case24464_goc",
"pglib_opf_case30000_goc",
"pglib_opf_case78484_epigrids",]=#

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



function generate_tex_opf(opf_results::Dict, coords; filename="benchmark_results_opf.tex")

    df_top = opf_results[:top]
    df_lifted_kkt = opf_results[:lifted_kkt]
    df_hybrid_kkt = opf_results[:hybrid_kkt]
    df_madncl = opf_results[:madncl]
    df_ma27 = opf_results[:ma27]
    df_ma86 = opf_results[:ma86]
    df_ma97 = opf_results[:ma97]
    

    methods = ["MadNLP+LiftedKKT (GPU)", "MadNLP+HybridKKT (GPU)", "MadNCL (GPU)",
            "Ipopt+Ma27 (CPU)","Ipopt+Ma86 (CPU)","Ipopt+Ma97 (CPU)"]
    subs = Dict(
        "MadNLP+LiftedKKT (GPU)" => [:iter, :soltime, :inittime, :adtime,
                          :lintime, :termination, :obj, :cvio],
        "MadNLP+HybridKKT (GPU)" => [:iter, :soltime, :inittime, :adtime,
                            :lintime, :termination, :obj, :cvio],
        "MadNCL (GPU)" => [:iter, :soltime, :inittime, :adtime,
                          :lintime, :termination, :obj, :cvio],
        "Ipopt+Ma27 (CPU)" => [:iter, :soltime, :adtime,
                          :termination, :obj, :cvio],
        "Ipopt+Ma86 (CPU)" => [:iter, :soltime, :adtime,
                          :termination, :obj, :cvio],
        "Ipopt+Ma97 (CPU)" => [:iter, :soltime, :adtime,
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
    raw_rows = Any[]
    for (i, row_top) in enumerate(eachrow(df_top))
        case = row_top.case_name
        clean_case = replace(case,
            r"^(api/|sad/)" => "",
            r"pglib_opf_case" => "",
            r"\.m$" => ""
        )
        row = Any[clean_case, format_k(row_top.nvar), format_k(row_top.ncon)]
        raw_row = Any[clean_case, row_top.nvar, row_top.ncon]

        methods = ["MadNLP+LiftedKKT (GPU)", "MadNLP+HybridKKT (GPU)", "MadNCL (GPU)",
            "Ipopt+Ma27 (CPU)","Ipopt+Ma86 (CPU)","Ipopt+Ma97 (CPU)"]
        for (df, method) in [(df_lifted_kkt, "MadNLP+LiftedKKT (GPU)"), (df_hybrid_kkt, "MadNLP+HybridKKT (GPU)"),
                            (df_madncl, "MadNCL (GPU)"), (df_ma27, "Ipopt+Ma27 (CPU)"),
                            (df_ma86, "Ipopt+Ma86 (CPU)"), (df_ma97, "Ipopt+Ma97 (CPU)")]
            df_row = df[i, :]
            for field in subs[method]
                val = get(df_row, field, missing)
                push!(row, format_val(field, val))
                push!(raw_row, val)
            end
        end

        push!(rows, row)
        push!(raw_rows, raw_row)
    end


    table_data = permutedims(reduce(hcat, rows))
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

    # Text file
    txt_filename = replace(filename, r"\.tex$" => ".txt")
    open(txt_filename, "w") do io
        pretty_table(
            io, table_data;
            header = (h_top, h_bottom),
            backend = Val(:text),
            alignment = :c
        )
    end

    # Raw CSV
    csv_filename = replace(filename, r"\.tex$" => ".csv")
    flat_header = vcat(["Case", "nvars", "ncons"], vcat([
        string(m, "_", f) for m in methods for f in subs[m]
    ]))
    df = DataFrame([Symbol(h) => col for (h, col) in zip(flat_header, eachcol(permutedims(reduce(hcat, raw_rows))))])
    CSV.write(csv_filename, df)


    # Full comparison
    selected = Dict(k => opf_results[k] for k in [:lifted_kkt, :hybrid_kkt, :madncl, :ma27, :ma86, :ma97])
    p = performance_profile(selected, df -> df.soltime)
    Plots.svg(p, replace(filename, r"\.tex$" => ""))

    small_list = Int64[]
    med_list= Int64[]
    large_list = Int64[]

    for (i, row) in enumerate(eachrow(opf_results[:top]))
        if row.nvar <= 2000
            push!(small_list, i)
        elseif row.nvar <= 20000
            push!(med_list, i)
        else
            push!(large_list, i)
        end
    end
    
    # Small: nvar < 2000
    ordered_keys = [:lifted_kkt, :hybrid_kkt, :madncl, :ma27, :ma86, :ma97]

    selected = OrderedDict(
        k => filter(row -> row.id in small_list, opf_results[k])
        for k in ordered_keys
        if haskey(opf_results, k)
    )
    # Now build an ordered list of Pairs
    if !isempty(selected[:lifted_kkt])
        p = performance_profile(selected, df -> df.soltime)
        Plots.svg(p, replace(filename, r"\.tex$" => "_small"))
    end


    # Medium: 2000 ≤ nvar ≤ 20000
    selected = OrderedDict(
        k => filter(row -> row.id in med_list, opf_results[k])
        for k in ordered_keys
        if haskey(opf_results, k)
    )
    # Now build an ordered list of Pairs
    if !isempty(selected[:lifted_kkt])
        p = performance_profile(selected, df -> df.soltime)
        Plots.svg(p, replace(filename, r"\.tex$" => "_small"))
    end

    # Large: nvar > 20000
    selected = OrderedDict(
        k => filter(row -> row.id in large_list, opf_results[k])
        for k in ordered_keys
        if haskey(opf_results, k)
    )
    # Now build an ordered list of Pairs
    if !isempty(selected[:lifted_kkt])
        p = performance_profile(selected, df -> df.soltime)
        Plots.svg(p, replace(filename, r"\.tex$" => "_small"))
    end


    # Log-log scatterplot of speedup vs nvar

    baseline = df_ma27
    n = nrow(df_top)

    scatter_data = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()

    for (method, df) in [
        ("MadNLP+LiftedKKT (GPU)", df_lifted_kkt),
        ("MadNLP+HybridKKT (GPU)", df_hybrid_kkt),
        ("MadNCL (GPU)", df_madncl),
        ("Ipopt+Ma86 (CPU)", df_ma86),
        ("Ipopt+Ma97 (CPU)", df_ma97)
    ]
        nv = Float64[]
        speedup = Float64[]
        for i in 1:n
            t_base = get(baseline[i, :], :soltime, missing)
            t = get(df[i, :], :soltime, missing)
            nv_i = get(df_top[i, :], :nvar, missing)

            if t !== missing && t_base !== missing && nv_i !== missing && t > 0 && t_base > 0
                push!(nv, float(nv_i))
                push!(speedup, t_base / t)
            end
        end
        scatter_data[method] = (nv, speedup)
    end

    p = plot(
        xlabel = "nvar", ylabel = "Speedup vs. Ma27",
        xscale = :log10, yscale = :log10,
        legend = :topleft, title = "Speedup vs. Problem Size",
        markerstrokewidth = 0
    )

    for (method, (nv, speedup)) in scatter_data
        scatter!(p, nv, speedup; label = method, ms=4)
    end

    svg_speedup = replace(filename, r"\.tex$" => "_speedup_vs_ma27.svg")
    savefig(p, svg_speedup)



end


function generate_tex_mpopf(mpopf_results, coords, curve_names; filename="benchmark_results_mpopf.tex")

    df_top = mpopf_results[:top]
    df_gpu_easy = mpopf_results[:gpu_easy]
    df_cpu_easy = mpopf_results[:cpu_easy]
    df_gpu_medium = mpopf_results[:gpu_medium]
    df_cpu_medium = mpopf_results[:cpu_medium]
    df_gpu_hard = mpopf_results[:gpu_hard]
    df_cpu_hard = mpopf_results[:cpu_hard]

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
    raw_rows = Any[]
    for (i, row_top) in enumerate(eachrow(df_top))
        case = row_top.case_name
        clean_case = replace(case, r"^pglib_opf_case" => "", r"\.m$" => "")

        row = Any[clean_case, format_k(row_top.nvar), format_k(row_top.ncon)]
        raw_row = Any[clean_case, row_top.nvar, row_top.ncon]

        #Easy
        row_gpu_easy = df_gpu_easy[i, :]
        for field in subs["GPU "*coords*" easy"]
            val = get(row_gpu_easy, field, missing)
            push!(row, format_val(field, val))
            push!(raw_row, val)
        end

        row_cpu_easy = df_cpu_easy[i, :]
        for field in subs["CPU "*coords*" easy"]
            val = get(row_cpu_easy, field, missing)
            push!(row, format_val(field, val))
            push!(raw_row, val)
        end

        #Medium
        row_gpu_medium = df_gpu_medium[i, :]
        for field in subs["GPU "*coords*" medium"]
            val = get(row_gpu_medium, field, missing)
            push!(row, format_val(field, val))
            push!(raw_row, val)
        end

        row_cpu_medium = df_cpu_medium[i, :]
        for field in subs["CPU "*coords*" medium"]
            val = get(row_cpu_medium, field, missing)
            push!(row, format_val(field, val))
            push!(raw_row, val)
        end

        #Hard
        row_gpu_hard = df_gpu_hard[i, :]
        for field in subs["GPU "*coords*" hard"]
            val = get(row_gpu_hard, field, missing)
            push!(row, format_val(field, val))
            push!(raw_row, val)
        end

        row_cpu_hard = df_cpu_hard[i, :]
        for field in subs["CPU "*coords*" hard"]
            val = get(row_cpu_hard, field, missing)
            push!(row, format_val(field, val))
            push!(raw_row, val)
        end


        push!(rows, row)
        push!(raw_rows, raw_row)
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


    # Write plain-text version (filename.txt)
    txt_filename = replace(filename, r"\.tex$" => ".txt")
    open(txt_filename, "w") do io
        pretty_table(
            io, table_data;
            header = (h_top, h_bottom),
            backend = Val(:text),
            alignment = :c
        )
    end

    # Write CSV version (raw values)
    csv_filename = replace(filename, r"\.tex$" => ".csv")
    flat_header = vcat(["Case", "nvars", "ncons"], vcat([
        string(m, "_", f) for m in methods for f in subs[m]
    ]))
    df = DataFrame([Symbol(h) => col for (h, col) in zip(flat_header, eachcol(permutedims(reduce(hcat, raw_rows))))])
    CSV.write(csv_filename, df)

    #Make charts
    selected = Dict(
        :gpu_easy => mpopf_results[:gpu_easy],
        :cpu_easy => mpopf_results[:cpu_easy],
    )
    p = performance_profile(selected, df -> df.soltime)
    Plots.svg(p, replace(filename, r"\.tex$" => "_easy"))

    selected = Dict(
        :gpu_medium => mpopf_results[:gpu_medium],
        :cpu_medium => mpopf_results[:cpu_medium],
    )
    p = performance_profile(selected, df -> df.soltime)
    Plots.svg(p, replace(filename, r"\.tex$" => "_medium"))

    selected = Dict(
        :gpu_hard => mpopf_results[:gpu_hard],
        :cpu_hard => mpopf_results[:cpu_hard],
    )
    p = performance_profile(selected, df -> df.soltime)
    Plots.svg(p, replace(filename, r"\.tex$" => "_hard"))


end

function generate_tex_stor_comp(stor_results, coords, comp_names; filename="benchmark_results_mpopf_storage.tex")

    df_top = stor_results[:top]
    df_gpu_no_cmp = stor_results[:gpu_no_cmp]
    df_gpu_cmp = stor_results[:gpu_cmp]
    df_gpu_nl_cmp = stor_results[:gpu_nl_cmp]
    df_cpu_no_cmp = stor_results[:cpu_no_cmp]
    df_cpu_cmp = stor_results[:cpu_cmp]
    df_cpu_nl_cmp = stor_results[:cpu_nl_cmp]

    # --- Build dynamic method names ---
    methods = String[]
    for comp in comp_names
        push!(methods, "GPU $coords $comp")
        push!(methods, "CPU $coords $comp")
    end

    # --- Define what fields each method has ---
    subs = Dict{String, Vector{Symbol}}()
    for method in methods
        if startswith(method, "GPU")
            subs[method] = [:nvar, :ncon, :iter, :soltime, :inittime, :adtime, :lintime, :termination, :obj, :cvio]
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
    raw_rows = Any[]
    for (i, row_top) in enumerate(eachrow(df_top))
        case = row_top.case_name
        clean_case = replace(case, r"^pglib_opf_case" => "", r"\.m$" => "")

        row = Any[clean_case]#, format_k(row_top.nvar), format_k(row_top.ncon)]
        raw_row = Any[clean_case]#, row_top.nvar, row_top.ncon]

        #Easy
        row_gpu_no_cmp = df_gpu_no_cmp[i, :]
        for field in subs["GPU "*coords*" no cc"]
            val = get(row_gpu_no_cmp, field, missing)
            if field in [:nvar, :ncon]
                push!(row, format_k(val))
            else
                push!(row, format_val(field, val))
            end
            push!(raw_row, val)
        end

        row_cpu_no_cmp = df_cpu_no_cmp[i, :]
        for field in subs["CPU "*coords*" no cc"]
            val = get(row_cpu_no_cmp, field, missing)
            push!(row, format_val(field, val))
            push!(raw_row, val)
        end

        #Medium
        row_gpu_cmp = df_gpu_cmp[i, :]
        for field in subs["GPU "*coords*" no cc"]
            val = get(row_gpu_cmp, field, missing)
            if field in [:nvar, :ncon]
                push!(row, format_k(val))
            else
                push!(row, format_val(field, val))
            end
            push!(raw_row, val)
        end

        row_cpu_cmp = df_cpu_cmp[i, :]
        for field in subs["CPU "*coords*" no cc"]
            val = get(row_cpu_cmp, field, missing)
            push!(row, format_val(field, val))
            push!(raw_row, val)
        end

        #Hard
        row_gpu_nl_cmp = df_gpu_nl_cmp[i, :]
        for field in subs["GPU "*coords*" no cc"]
            val = get(row_gpu_nl_cmp, field, missing)
            if field in [:nvar, :ncon]
                push!(row, format_k(val))
            else
                push!(row, format_val(field, val))
            end
            push!(raw_row, val)
        end

        row_cpu_nl_cmp = df_cpu_nl_cmp[i, :]
        for field in subs["CPU "*coords*" no cc"]
            val = get(row_cpu_nl_cmp, field, missing)
            push!(row, format_val(field, val))
            push!(raw_row, val)
        end


        push!(rows, row)
        push!(raw_rows, raw_row)
    end

    table_data = permutedims(reduce(hcat, rows))

    # --- Header construction ---
    h_top    = ["Case"]
    h_bottom = [""]

    for m in methods
        n = length(subs[m])
        push!(h_top, m)
        append!(h_top, fill("", n-1))
        append!(h_bottom, string.(subs[m]))
    end

    # --- Group boundary vlines ---
    function group_boundaries(methods, subs)
        idx = Int[0, 1]  # vertical lines at case label and after header
        col = 1
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


    # Write plain-text version (filename.txt)
    txt_filename = replace(filename, r"\.tex$" => ".txt")
    open(txt_filename, "w") do io
        pretty_table(
            io, table_data;
            header = (h_top, h_bottom),
            backend = Val(:text),
            alignment = :c
        )
    end

    # Write CSV version (raw values)
    csv_filename = replace(filename, r"\.tex$" => ".csv")
    flat_header = vcat(["Case"], vcat([
        string(m, "_", f) for m in methods for f in subs[m]
    ]))
    df = DataFrame([Symbol(h) => col for (h, col) in zip(flat_header, eachcol(permutedims(reduce(hcat, raw_rows))))])
    CSV.write(csv_filename, df)

    #Make charts
    selected = Dict(
        :gpu_no_cmp => stor_results[:gpu_no_cmp],
        :cpu_no_cmp => stor_results[:cpu_no_cmp],
    )
    p = performance_profile(selected, df -> df.soltime)
    Plots.svg(p, replace(filename, r"\.tex$" => "_no_cmp"))

    selected = Dict(
        :gpu_cmp => stor_results[:gpu_cmp],
        :cpu_cmp => stor_results[:cpu_cmp],
    )
    p = performance_profile(selected, df -> df.soltime)
    Plots.svg(p, replace(filename, r"\.tex$" => "_cmp"))

    selected = Dict(
        :gpu_nl_cmp => stor_results[:gpu_nl_cmp],
        :cpu_nl_cmp => stor_results[:cpu_nl_cmp],
    )
    p = performance_profile(selected, df -> df.soltime)
    Plots.svg(p, replace(filename, r"\.tex$" => "_nl_cmp"))


end


function solve_static_cases(cases, tol, coords; case_style = "default")

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
    ~ = MadNCL.madncl(model_gpu, tol = tol, max_iter = 3, disable_garbage_collector=true, dual_initialized=true)
    ~ = madnlp(model_gpu, tol = tol, max_iter = 3, disable_garbage_collector=true, dual_initialized=true, linear_solver=MadNLPGPU.CUDSSSolver,
                cudss_algorithm=MadNLP.LDL,
                kkt_system=HybridKKT.HybridCondensedKKTSystem,
                equality_treatment=MadNLP.EnforceEquality,
                fixed_variable_treatment=MadNLP.MakeParameter,)

    model_cpu, ~ = opf_model("pglib_opf_case3_lmbd"; form=form)
    ~ = ipopt(model_cpu, tol = tol, max_iter = 3, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma27")
    ~ = ipopt(model_cpu, tol = tol, max_iter = 3, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma86")
    ~ = ipopt(model_cpu, tol = tol, max_iter = 3, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma97")



    df_top = DataFrame(nvar = Int[], ncon = Int[], case_name = String[])

    df_lifted_kkt = DataFrame(id = Int[], iter = Float64[], soltime = Float64[], inittime = Float64[], adtime = Float64[], lintime = Float64[], 
        termination = String[], obj = Float64[], cvio = Float64[])
    df_hybrid_kkt = similar(df_lifted_kkt)
    df_madncl = similar(df_lifted_kkt)

    df_ma27 = DataFrame(id = Int[], iter = Int[], soltime = Float64[], adtime = Float64[], termination = String[], obj = Float64[], cvio = Float64[])
    df_ma86 = similar(df_ma27)
    df_ma97 = similar(df_ma27)

    for (i, case) in enumerate(cases)
        println(case)

        if case_style == "default"
            case = case*".m"
        elseif case_style == "api"
            case = "api/"*case*"__api.m"
        elseif case_style == "sad"
            case = "sad/"*case*"__sad.m"
        else
            error("Invalid case style")
        end


        #GPU 
        m_gpu, v_gpu, c_gpu = opf_model(case; backend = CUDABackend(), form=form)   

        result_lifted_kkt = madnlp(m_gpu, tol=tol, max_wall_time = max_wall_time, disable_garbage_collector=true, dual_initialized=true)
        c = evaluate(m_gpu, result_lifted_kkt)
        push!(df_lifted_kkt, (i, result_lifted_kkt.counters.k, result_lifted_kkt.counters.total_time, result_lifted_kkt.counters.init_time, result_lifted_kkt.counters.eval_function_time, 
        result_lifted_kkt.counters.linear_solver_time, termination_code(result_lifted_kkt.status), result_lifted_kkt.objective, c))
        push!(df_top, (m_gpu.meta.nvar, m_gpu.meta.ncon, case))

        result_hybrid_kkt = madnlp(m_gpu, tol=tol, max_wall_time = max_wall_time, disable_garbage_collector=true, dual_initialized=true, linear_solver=MadNLPGPU.CUDSSSolver,
                                    cudss_algorithm=MadNLP.LDL,
                                    kkt_system=HybridKKT.HybridCondensedKKTSystem,
                                    equality_treatment=MadNLP.EnforceEquality,
                                    fixed_variable_treatment=MadNLP.MakeParameter,)
        c = evaluate(m_gpu, result_hybrid_kkt)
        push!(df_hybrid_kkt, (i, result_hybrid_kkt.counters.k, result_hybrid_kkt.counters.total_time, result_hybrid_kkt.counters.init_time, result_hybrid_kkt.counters.eval_function_time, 
        result_hybrid_kkt.counters.linear_solver_time, termination_code(result_hybrid_kkt.status), result_hybrid_kkt.objective, c))

        result_madncl = MadNCL.madncl(m_gpu, tol=tol, max_wall_time = max_wall_time, disable_garbage_collector=true, dual_initialized=true)
        c = evaluate(m_gpu, result_madncl)
        push!(df_madncl, (i, result_madncl.counters.k, result_madncl.counters.total_time, result_madncl.counters.init_time, result_madncl.counters.eval_function_time, 
        result_madncl.counters.linear_solver_time, termination_code(result_madncl.status), result_madncl.objective, c))


        #CPU
        m_cpu, v_cpu, c_cpu = opf_model(case; form=form)

        result_ma27 = ipopt(m_cpu, tol = tol, max_wall_time=max_wall_time, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma27", honor_original_bounds = "no", print_timing_statistics = "yes", output_file = "ipopt_output")
        it, tot, ad = ipopt_stats("ipopt_output")
        c = evaluate(m_cpu, result_ma27)
        push!(df_ma27, (i, it, tot, ad, termination_code(result_ma27.solver_specific[:internal_msg]), result_ma27.objective, c))

        result_ma86 = ipopt(m_cpu, tol = tol, max_wall_time=max_wall_time, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma86", honor_original_bounds = "no", print_timing_statistics = "yes", output_file = "ipopt_output")
        it, tot, ad = ipopt_stats("ipopt_output")
        c = evaluate(m_cpu, result_ma86)
        push!(df_ma86, (i, it, tot, ad, termination_code(result_ma86.solver_specific[:internal_msg]), result_ma86.objective, c))

        result_ma97 = ipopt(m_cpu, tol = tol, max_wall_time=max_wall_time, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma97", honor_original_bounds = "no", print_timing_statistics = "yes", output_file = "ipopt_output")
        it, tot, ad = ipopt_stats("ipopt_output")
        c = evaluate(m_cpu, result_ma97)
        push!(df_ma97, (i, it, tot, ad, termination_code(result_ma97.solver_specific[:internal_msg]), result_ma97.objective, c))
    end

    opf_results = Dict(:top => df_top,
                :lifted_kkt => df_lifted_kkt,
                :hybrid_kkt => df_hybrid_kkt,
                :madncl => df_madncl,
                :ma27 => df_ma27,
                :ma86 => df_ma86,
                :ma97 => df_ma97)
    
    generate_tex_opf(opf_results, coords; filename = "benchmark_results_opf_" * case_style * "_tol_" * replace(@sprintf("%.0e", 1 / tol), r"\+0?" => "")*".tex")

    return opf_results
end

curves = Dict("easy" => [.64, .60, .58, .56, .56, .58, .64, .76, .87, .95, .99, 1.0, .99, 1.0, 1.0,
    .97, .96, .96, .93, .92, .92, .93, .87, .72, .64],
    "medium" => [.88, .90, .88, .86, .87, .88, .9, .92, .93, .95, .97, 1.0, .99, 1.0, 1.0,
    .97, .96, .96, .93, .92, .92, .93, .89, .85, .82],
    "hard" => [.52, .60, .53, .59, .51, .62, .65, .76, .87, .95, .99, 1.01, .99, 1.0, 1.02,
    .92, 1.0, .9, .93, .84, .92, .93, .85, .73, .62])


function solve_mp_cases(cases, curves, tol, coords; case_style = "default", storage = false)

    

    if coords == "Polar"
        form = :polar
    elseif coords == "Rectangular"
        form = :rect
    else
        error("Wrong coords")
    end

    df_top = DataFrame(nvar = Int[], ncon = Int[], case_name = String[])

    df_lifted_kkt_easy = DataFrame(id = Int[], iter = Float64[], soltime = Float64[], inittime = Float64[], adtime = Float64[], lintime = Float64[], 
        termination = String[], obj = Float64[], cvio = Float64[])

    df_lifted_kkt_medium = similar(df_lifted_kkt_easy)
    df_lifted_kkt_hard = similar(df_lifted_kkt_easy)

    df_hybrid_kkt_easy = similar(df_lifted_kkt_easy)
    df_hybrid_kkt_medium = similar(df_lifted_kkt_easy)
    df_hybrid_kkt_hard = similar(df_lifted_kkt_easy)

    df_madncl_easy = similar(df_lifted_kkt_easy)
    df_madncl_medium = similar(df_lifted_kkt_easy)
    df_madncl_hard = similar(df_lifted_kkt_easy)

    df_ma27_easy = DataFrame(id = Int[], iter = Int[], soltime = Float64[], adtime = Float64[], termination = String[], obj = Float64[], cvio = Float64[])
    df_ma27_medium = similar(df_ma27_easy)
    df_ma27_hard = similar(df_ma27_easy)

    df_ma86_easy = similar(df_ma27_easy)
    df_ma86_medium = similar(df_ma27_easy)
    df_ma86_hard = similar(df_ma27_easy)

    df_ma97_easy = similar(df_ma27_easy)
    df_ma97_medium = similar(df_ma27_easy)
    df_ma97_hard = similar(df_ma27_easy)


    #Compile time on smallest case
    if !storage
        test_case = "pglib_opf_case3_lmbd"
        max_wall_time = Float64(2500)
    else
        test_case = "pglib_opf_case3_lmbd_storage.m"
        max_wall_time = Float64(4500)
    end
    model_gpu, ~ = mpopf_model(test_case, [1,1,1]; backend = CUDABackend(), form=form)
    ~ = madnlp(model_gpu, tol = tol, max_iter = 3, disable_garbage_collector=true, dual_initialized=true)
    ~ = MadNCL.madncl(model_gpu, tol = tol, max_iter = 3, disable_garbage_collector=true, dual_initialized=true)
    ~ = madnlp(model_gpu, tol = tol, max_iter = 3, disable_garbage_collector=true, dual_initialized=true, linear_solver=MadNLPGPU.CUDSSSolver,
                cudss_algorithm=MadNLP.LDL,
                kkt_system=HybridKKT.HybridCondensedKKTSystem,
                equality_treatment=MadNLP.EnforceEquality,
                fixed_variable_treatment=MadNLP.MakeParameter,)

    model_cpu, ~ = mpopf_model(test_case, [1,1,1]; form=form)
    ~ = ipopt(model_cpu, tol = tol, max_iter = 3, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma27")
    ~ = ipopt(model_cpu, tol = tol, max_iter = 3, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma86")
    ~ = ipopt(model_cpu, tol = tol, max_iter = 3, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma97")


    for (i, case) in enumerate(cases)
        if !storage
            if case_style == "default"
                case = case*".m"
            elseif case_style == "api"
                case = "api/"*case*"__api.m"
            elseif case_style == "sad"
                case = "sad/"*case*"__sad.m"
            else
                error("Invalid case style")
            end
        else
            if case_style == "default"
                case = case*"_storage.m"
            elseif case_style == "api"
                case = case*"__api_storage.m"
            elseif case_style == "sad"
                case = case*"__sad_storage.m"
            else
                error("Invalid case style")
            end
        end
        for (curve_name, curve) in curves            

            #GPU
            m_gpu, v_gpu, c_gpu = mpopf_model(case, curve; backend = CUDABackend(), form=form)   
            result_lifted_kkt = madnlp(m_gpu, tol=tol, max_wall_time = max_wall_time, disable_garbage_collector=true, dual_initialized=true)

            c = evaluate(m_gpu, result_lifted_kkt)        

            row_lifted_kkt = (i, result_lifted_kkt.counters.k, result_lifted_kkt.counters.total_time, result_lifted_kkt.counters.init_time, result_lifted_kkt.counters.eval_function_time,
                result_lifted_kkt.counters.linear_solver_time, termination_code(result_lifted_kkt.status), result_lifted_kkt.objective, c)

            result_hybrid_kkt = madnlp(m_gpu, tol=tol, max_wall_time = max_wall_time, disable_garbage_collector=true, dual_initialized=true, linear_solver=MadNLPGPU.CUDSSSolver,
                                    cudss_algorithm=MadNLP.LDL,
                                    kkt_system=HybridKKT.HybridCondensedKKTSystem,
                                    equality_treatment=MadNLP.EnforceEquality,
                                    fixed_variable_treatment=MadNLP.MakeParameter,)
            c = evaluate(m_gpu, result_hybrid_kkt)
            row_hybrid_kkt = (i, result_hybrid_kkt.counters.k, result_hybrid_kkt.counters.total_time, result_hybrid_kkt.counters.init_time, result_hybrid_kkt.counters.eval_function_time,
                result_hybrid_kkt.counters.linear_solver_time, termination_code(result_hybrid_kkt.status), result_hybrid_kkt.objective, c)

            result_madncl = MadNCL.madncl(m_gpu, tol=tol, max_wall_time = max_wall_time, disable_garbage_collector=true, dual_initialized=true)
            c = evaluate(m_gpu, result_madncl)
            row_madncl = (i, result_madncl.counters.k, result_madncl.counters.total_time, result_madncl.counters.init_time, result_madncl.counters.eval_function_time,
                result_madncl.counters.linear_solver_time, termination_code(result_madncl.status), result_madncl.objective, c)



            #CPU
            m_cpu, v_cpu, c_cpu = mpopf_model(case, curve; form = form)

            result_ma27 = ipopt(m_cpu, tol = tol, max_wall_time=max_wall_time, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma27", honor_original_bounds = "no", print_timing_statistics = "yes", output_file = "ipopt_output")
            it, tot, ad = ipopt_stats("ipopt_output")
            c = evaluate(m_cpu, result_ma27)
            row_ma27 = (i, it, tot, ad, termination_code(result_ma27.solver_specific[:internal_msg]),  result_ma27.objective, c)

            result_ma86 = ipopt(m_cpu, tol = tol, max_wall_time=max_wall_time, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma86", honor_original_bounds = "no", print_timing_statistics = "yes", output_file = "ipopt_output")
            it, tot, ad = ipopt_stats("ipopt_output")
            c = evaluate(m_cpu, result_ma86)
            row_ma86 = (i, it, tot, ad, termination_code(result_ma86.solver_specific[:internal_msg]),  result_ma86.objective, c)

            result_ma97 = ipopt(m_cpu, tol = tol, max_wall_time=max_wall_time, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma97", honor_original_bounds = "no", print_timing_statistics = "yes", output_file = "ipopt_output")
            it, tot, ad = ipopt_stats("ipopt_output")
            c = evaluate(m_cpu, result_ma97)
            row_ma97 = (i, it, tot, ad, termination_code(result_ma97.solver_specific[:internal_msg]),  result_ma97.objective, c)
            
            if curve_name == "easy"
                push!(df_top, (m_gpu.meta.nvar, m_gpu.meta.ncon, case))
                push!(df_lifted_kkt_easy, row_lifted_kkt)
                push!(df_hybrid_kkt_easy, row_hybrid_kkt)
                push!(df_madncl_easy, row_madncl)
                push!(df_ma27_easy, row_ma27)
                push!(df_ma86_easy, row_ma86)
                push!(df_ma97_easy, row_ma97)
            elseif curve_name == "medium"
                push!(df_lifted_kkt_medium, row_lifted_kkt)
                push!(df_hybrid_kkt_medium, row_hybrid_kkt)
                push!(df_madncl_medium, row_madncl)
                push!(df_ma27_medium, row_ma27)
                push!(df_ma86_medium, row_ma86)
                push!(df_ma97_medium, row_ma97)
            elseif curve_name == "hard"
                push!(df_lifted_kkt_hard, row_lifted_kkt)
                push!(df_hybrid_kkt_hard, row_hybrid_kkt)
                push!(df_madncl_hard, row_madncl)
                push!(df_ma27_hard, row_ma27)
                push!(df_ma86_hard, row_ma86)
                push!(df_ma97_hard, row_ma97)
            end

            
        end

    end

    curve_names = collect(keys(curves))

    mpopf_results = Dict(:top => df_top,
                :lifted_kkt_easy => df_lifted_kkt_easy,
                :hybrid_kkt_easy => df_hybrid_kkt_easy,
                :madncl_easy => df_madncl_easy,
                :lifted_kkt_medium => df_lifted_kkt_medium,
                :hybrid_kkt_medium => df_hybrid_kkt_medium,
                :madncl_medium => df_madncl_medium,
                :lifted_kkt_hard => df_lifted_kkt_hard,
                :hybrid_kkt_hard => df_hybrid_kkt_hard,
                :madncl_hard => df_madncl_hard,
                :ma27_easy => df_ma27_easy,
                :ma86_easy => df_ma86_easy,
                :ma97_easy => df_ma97_easy,
                :ma27_medium => df_ma27_medium,
                :ma86_medium => df_ma86_medium,
                :ma97_medium => df_ma97_medium,
                :ma27_hard => df_ma27_hard,
                :ma86_hard => df_ma86_hard,
                :ma97_hard => df_ma97_hard,)
    
    #if !storage
    #    generate_tex_mpopf(mpopf_results, coords, curve_names; filename="benchmark_results_mpopf_" * case_style*"_tol_" * replace(@sprintf("%.0e", 1 / tol), r"\+0?" => "")*".tex")
    #else
    #    generate_tex_mpopf(mpopf_results, coords, curve_names; filename="benchmark_results_mpopf_storage_" * case_style*"_tol_" * replace(@sprintf("%.0e", 1 / tol), r"\+0?" => "")*".tex")
    #end
    return mpopf_results
end


function solve_stor_cases_comp(cases, tol, coords; case_style = "default", curve = [.88, .90, .88, .86, .87, .88, .9, .92, .93, .95, .97, 1.0, .99, 1.0, 1.0,
    .97, .96, .96, .93, .92, .92, .93, .89, .85, .82])

    max_wall_time = Float64(4500)

    function example_func(d, srating)
        return d + .2/srating*d^2
    end

    if coords == "Polar"
        form = :polar
    elseif coords == "Rectangular"
        form = :rect
    else
        error("Wrong coords")
    end

    df_top = DataFrame(case_name = String[])

    df_gpu_no_cmp = DataFrame(nvar = Int[], ncon = Int[], id = Int[], iter = Float64[], soltime = Float64[], inittime = Float64[], adtime = Float64[], lintime = Float64[], 
        termination = String[], obj = Float64[], cvio = Float64[])
    df_gpu_cmp = similar(df_gpu_no_cmp)
    df_gpu_nl_cmp = similar(df_gpu_no_cmp)

    df_cpu_no_cmp = DataFrame(id = Int[], iter = Int[], soltime = Float64[], adtime = Float64[], termination = String[], obj = Float64[], cvio = Float64[])
    df_cpu_cmp = similar(df_cpu_no_cmp)
    df_cpu_nl_cmp = similar(df_cpu_no_cmp)


    #Compile time on smallest case
    model_gpu, ~ = mpopf_model("pglib_opf_case3_lmbd_storage.m", [1,1,1]; backend = CUDABackend(), form=form)
    ~ = madnlp(model_gpu, tol = tol, max_iter = 3, disable_garbage_collector=true, dual_initialized=true)
    model_gpu, ~ = mpopf_model("pglib_opf_case3_lmbd_storage.m", [1,1,1]; backend = CUDABackend(), form=form, storage_complementarity_constraint = true)
    ~ = madnlp(model_gpu, tol = tol, max_iter = 3, disable_garbage_collector=true, dual_initialized=true)
    model_gpu, ~ = mpopf_model("pglib_opf_case3_lmbd_storage.m", [1,1,1], example_func; backend = CUDABackend(), form=form)
    ~ = madnlp(model_gpu, tol = tol, max_iter = 3, disable_garbage_collector=true, dual_initialized=true)

    model_cpu, ~ = mpopf_model("pglib_opf_case3_lmbd_storage.m", [1,1,1]; form=form)
    ~ = ipopt(model_cpu, tol = tol, max_iter = 3, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma27")
    model_cpu, ~ = mpopf_model("pglib_opf_case3_lmbd_storage.m", [1,1,1]; form=form, storage_complementarity_constraint = true)
    ~ = ipopt(model_cpu, tol = tol, max_iter = 3, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma27")
    model_cpu, ~ = mpopf_model("pglib_opf_case3_lmbd_storage.m", [1,1,1], example_func; form=form)
    ~ = ipopt(model_cpu, tol = tol, max_iter = 3, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma27")


    for (i, case) in enumerate(cases)

        #Requires modified to be created, with "_storage.m" appended
        if case_style == "default"
            case = case*"_storage.m"
        elseif case_style == "api"
            case = case*"__api_storage.m"
        elseif case_style == "sad"
            case = case*"__sad_storage.m"
        else
            error("Invalid case style")
        end
                 

        
        #No complementarity constraint
        m_gpu, v_gpu, c_gpu = mpopf_model(case, curve; backend = CUDABackend(), form=form)   
        result_gpu = madnlp(m_gpu, tol=tol, max_wall_time = max_wall_time, disable_garbage_collector=true, dual_initialized=true)

        c = evaluate(m_gpu, result_gpu)        

        row_gpu = (m_gpu.meta.nvar, m_gpu.meta.ncon, i, result_gpu.counters.k, result_gpu.counters.total_time, result_gpu.counters.init_time, result_gpu.counters.eval_function_time,
            result_gpu.counters.linear_solver_time, termination_code(result_gpu.status), result_gpu.objective, c)
        
        push!(df_top, (case,))
        push!(df_gpu_no_cmp, row_gpu)

        m_cpu, v_cpu, c_cpu = mpopf_model(case, curve; form = form)
        result_cpu = ipopt(m_cpu, tol = tol, max_wall_time=max_wall_time, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma27", honor_original_bounds = "no", print_timing_statistics = "yes", output_file = "ipopt_output")

        it, tot, ad = ipopt_stats("ipopt_output")
        c = evaluate(m_cpu, result_cpu)
        row_cpu = (i, it, tot, ad, termination_code(result_cpu.solver_specific[:internal_msg]),  result_cpu.objective, c)
        push!(df_cpu_no_cmp, row_cpu)

        #Complementarity constraint enforced
        m_gpu, v_gpu, c_gpu = mpopf_model(case, curve; backend = CUDABackend(), form=form, storage_complementarity_constraint = true)   
        result_gpu = madnlp(m_gpu, tol=tol, max_wall_time = max_wall_time, disable_garbage_collector=true, dual_initialized=true)

        c = evaluate(m_gpu, result_gpu)        

        row_gpu = (m_gpu.meta.nvar, m_gpu.meta.ncon, i, result_gpu.counters.k, result_gpu.counters.total_time, result_gpu.counters.init_time, result_gpu.counters.eval_function_time,
            result_gpu.counters.linear_solver_time, termination_code(result_gpu.status), result_gpu.objective, c)
        
        push!(df_gpu_cmp, row_gpu)

        m_cpu, v_cpu, c_cpu = mpopf_model(case, curve; form = form, storage_complementarity_constraint=true)
        result_cpu = ipopt(m_cpu, tol = tol, max_wall_time=max_wall_time, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma27", honor_original_bounds = "no", print_timing_statistics = "yes", output_file = "ipopt_output")

        it, tot, ad = ipopt_stats("ipopt_output")
        c = evaluate(m_cpu, result_cpu)
        row_cpu = (i, it, tot, ad, termination_code(result_cpu.solver_specific[:internal_msg]),  result_cpu.objective, c)
        push!(df_cpu_cmp, row_cpu)

        #Replace piecewise complementarity constraint with NL smooth function
        m_gpu, v_gpu, c_gpu = mpopf_model(case, curve, example_func; backend = CUDABackend(), form=form)   
        result_gpu = madnlp(m_gpu, tol=tol, max_wall_time = max_wall_time, disable_garbage_collector=true, dual_initialized=true)

        c = evaluate(m_gpu, result_gpu)        

        row_gpu = (m_gpu.meta.nvar, m_gpu.meta.ncon, i, result_gpu.counters.k, result_gpu.counters.total_time, result_gpu.counters.init_time, result_gpu.counters.eval_function_time,
            result_gpu.counters.linear_solver_time, termination_code(result_gpu.status), result_gpu.objective, c)
        
        push!(df_gpu_nl_cmp, row_gpu)

        m_cpu, v_cpu, c_cpu = mpopf_model(case, curve, example_func; form = form)
        result_cpu = ipopt(m_cpu, tol = tol, max_wall_time=max_wall_time, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma27", honor_original_bounds = "no", print_timing_statistics = "yes", output_file = "ipopt_output")

        it, tot, ad = ipopt_stats("ipopt_output")
        c = evaluate(m_cpu, result_cpu)
        row_cpu = (i, it, tot, ad, termination_code(result_cpu.solver_specific[:internal_msg]),  result_cpu.objective, c)
        push!(df_cpu_nl_cmp, row_cpu)
           
    end

    comp_names = ["no cc", "cc", "nl cc"]

    stor_results = Dict(:top => df_top,
                :gpu_no_cmp => df_gpu_no_cmp,
                :gpu_cmp => df_gpu_cmp,
                :gpu_nl_cmp => df_gpu_nl_cmp,
                :cpu_no_cmp => df_cpu_no_cmp,
                :cpu_cmp => df_cpu_cmp,
                :cpu_nl_cmp => df_cpu_nl_cmp,)
    
    generate_tex_stor_comp(stor_results, coords, comp_names; filename="benchmark_results_mpopf_storage_comps_" * case_style*"_tol_" * replace(@sprintf("%.0e", 1 / tol), r"\+0?" => "")*".tex")
    return stor_results
end

sc_cases = [("data/C3E4N00073D1_scenario_303.json", "data/C3E4N00073D1_scenario_303_solution.json"),
            ("data/C3E4N00073D3_scenario_303.json", "data/C3E4N00073D3_scenario_303_solution.json")]

function solve_sc_cases(cases, tol, include_ctg)

    df_top = DataFrame(nvar = Int[], ncon = Int[], case_name = String[])

    df_gpu = DataFrame(id = Int[], iter = Float64[], soltime = Float64[], inittime = Float64[], adtime = Float64[], lintime = Float64[], 
        termination = String[], obj = Float64[], cvio = Float64[])

    df_cpu = DataFrame(id = Int[], iter = Int[], soltime = Float64[], adtime = Float64[], termination = String[], obj = Float64[], cvio = Float64[])


    #Compile time on smallest case
    
    test_case = "data/C3E4N00073D1_scenario_303.json"
    test_uc_case = "data/C3E4N00073D1_scenario_303_solution.json"
    max_wall_time = Float64(7000)
    
    model_gpu, ~ = scopf_model(test_case, test_uc_case; backend = CUDABackend(), include_ctg = include_ctg)
    ~ = madnlp(model_gpu, tol = tol, max_iter = 3, disable_garbage_collector=true, dual_initialized=true)

    model_cpu, ~ = scopf_model(test_case, test_uc_case; include_ctg = include_ctg)
    ~ = ipopt(model_cpu, tol = tol, max_iter = 3, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma27")


    for (i, case) in enumerate(cases)          

        #GPU
        (problem_case, uc_case) = case
        m_gpu, v_gpu, c_gpu = scopf_model(problem_case, uc_case; backend = CUDABackend(), include_ctg = include_ctg)   
        result_gpu = madnlp(m_gpu, tol=tol, max_wall_time = max_wall_time, disable_garbage_collector=true, dual_initialized=true)

        c = evaluate(m_gpu, result_gpu)        

        row_gpu = (i, result_gpu.counters.k, result_gpu.counters.total_time, result_gpu.counters.init_time, result_gpu.counters.eval_function_time,
            result_gpu.counters.linear_solver_time, termination_code(result_gpu.status), result_gpu.objective, c)

        m_cpu, v_cpu, c_cpu = scopf_model(problem_case, uc_case); include_ctg = include_ctg
        result_cpu = ipopt(m_cpu, tol = tol, max_wall_time=max_wall_time, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = tol, linear_solver = "ma27", honor_original_bounds = "no", print_timing_statistics = "yes", output_file = "ipopt_output")

        it, tot, ad = ipopt_stats("ipopt_output")

        c = evaluate(m_cpu, result_cpu)
        
        row_cpu = (i, it, tot, ad, termination_code(result_cpu.solver_specific[:internal_msg]),  result_cpu.objective, c)

        push!(df_top, (m_gpu.meta.nvar, m_gpu.meta.ncon, case))
        push!(df_gpu, row_gpu)
        push!(df_cpu, row_cpu)
    end

    scopf_results = Dict(:top => df_top,
                :gpu => df_gpu
                :cpu => df_cpu)
 
    #generate_tex_mpopf(mpopf_results, coords, curve_names; filename="benchmark_results_mpopf_" * case_style*"_tol_" * replace(@sprintf("%.0e", 1 / tol), r"\+0?" => "")*".tex")
   
    return scopf_results
end