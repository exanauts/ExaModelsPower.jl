convert_data(data::N, backend) where {names,N<:NamedTuple{names}} =
    NamedTuple{names}(convert_array(d, backend) for d in data)

function parse_ac_power_data(filename)
    d, f = splitdir(filename)
    name, ext = splitext(f)

    if isfile(joinpath(TMPDIR, name) * ".jld2")
        @info "Loading cached JLD2 file"
        loaded = JLD2.load(joinpath(TMPDIR, name) * ".jld2")
        return loaded["data"], loaded["dicts"]
    else
        ff = if isfile(filename)
            filename
        elseif isfile(joinpath(TMPDIR, name) * ".m")
            joinpath(TMPDIR, name) * ".m"
        else
            joinpath(PGLib.PGLib_opf, name * ".m")
        end
        @info "Loading MATPOWER file"
        return process_ac_power_data(ff)
    end
end

function process_ac_power_data(filename)
    data = PowerModels.parse_file(filename)
    PowerModels.standardize_cost_terms!(data, order = 2)
    PowerModels.calc_thermal_limits!(data)

    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    dicts = (
        arc = Dict(a => k for (k, a) in enumerate(ref[:arcs])),
        bus = Dict(k => i for (i, (k, v)) in enumerate(ref[:bus])),
        gen = Dict(k => i for (i, (k, v)) in enumerate(ref[:gen])),
        branch = Dict(k => i for (i, (k, v)) in enumerate(ref[:branch])),
    )

    data =  (
        baseMVA = [ref[:baseMVA]],
        bus = [
            begin
                bus_loads = [ref[:load][l] for l in ref[:bus_loads][k]]
                bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][k]]
                pd = sum(load["pd"] for load in bus_loads; init = 0.0)
                gs = sum(shunt["gs"] for shunt in bus_shunts; init = 0.0)
                qd = sum(load["qd"] for load in bus_loads; init = 0.0)
                bs = sum(shunt["bs"] for shunt in bus_shunts; init = 0.0)
                (i = dicts.bus[k], pd = pd, gs = gs, qd = qd, bs = bs, bus_type = v["bus_type"])
            end for (k, v) in ref[:bus]
                ],
        gen = [
            (
                i = dicts.gen[k],
                cost1 = v["cost"][1],
                cost2 = v["cost"][2],
                cost3 = v["cost"][3],
                bus = dicts.bus[v["gen_bus"]],
            ) for (k, v) in ref[:gen]
        ],
        arc = [
            (i = k, rate_a = ref[:branch][l]["rate_a"], bus = dicts.bus[i]) for
            (k, (l, i, j)) in enumerate(ref[:arcs])
                ],
        branch = [
            begin
                f_idx = dicts.arc[i, branch["f_bus"], branch["t_bus"]]
                t_idx = dicts.arc[i, branch["t_bus"], branch["f_bus"]]
                g, b = PowerModels.calc_branch_y(branch)
                tr, ti = PowerModels.calc_branch_t(branch)
                ttm = tr^2 + ti^2
                g_fr = branch["g_fr"]
                b_fr = branch["b_fr"]
                g_to = branch["g_to"]
                b_to = branch["b_to"]
                c1 = (-g * tr - b * ti) / ttm
                c2 = (-b * tr + g * ti) / ttm
                c3 = (-g * tr + b * ti) / ttm
                c4 = (-b * tr - g * ti) / ttm
                c5 = (g + g_fr) / ttm
                c6 = (b + b_fr) / ttm
                c7 = (g + g_to)
                c8 = (b + b_to)
                (
                    i = dicts.branch[i],
                    j = 1,
                    f_idx = f_idx,
                    t_idx = t_idx,
                    f_bus = dicts.bus[branch["f_bus"]],
                    t_bus = dicts.bus[branch["t_bus"]],
                    c1 = c1,
                    c2 = c2,
                    c3 = c3,
                    c4 = c4,
                    c5 = c5,
                    c6 = c6,
                    c7 = c7,
                    c8 = c8,
                    rate_a_sq = branch["rate_a"]^2,
                )
            end for (i, branch) in ref[:branch]
                ],
        storage = isempty(ref[:storage]) ?  empty_data = Vector{NamedTuple{(:i,), Tuple{Int64}}}() : [
            begin
                (c = i,
                Einit = stor["energy"],
                etac = stor["charge_efficiency"],
                etad = stor["discharge_efficiency"],
                Srating = stor["thermal_rating"],
                Zr = stor["r"],
                Zim = stor["x"],
                Pexts = stor["ps"],
                Qexts = stor["qs"],
                bus = dicts.bus[stor["storage_bus"]])
            end for (i, stor) in ref[:storage]
        ],
        ref_buses = [dicts.bus[i] for (i, k) in ref[:ref_buses]],
        vmax = [v["vmax"] for (k, v) in ref[:bus]],
        vmin = [v["vmin"] for (k, v) in ref[:bus]],
        pmax = [v["pmax"] for (k, v) in ref[:gen]],
        pmin = [v["pmin"] for (k, v) in ref[:gen]],
        qmax = [v["qmax"] for (k, v) in ref[:gen]],
        qmin = [v["qmin"] for (k, v) in ref[:gen]],
        rate_a = [ref[:branch][l]["rate_a"] for (k, (l, i, j)) in enumerate(ref[:arcs])],
        angmax = [b["angmax"] for (i, b) in ref[:branch]],
        angmin = [b["angmin"] for (i, b) in ref[:branch]],
        vm0 = [v["vm"] for (k, v) in ref[:bus]],
        va0 = [v["va"] for (k, v) in ref[:bus]],
        pg0 = [v["pg"] for (k, v) in ref[:gen]],
        qg0 = [v["qg"] for (k, v) in ref[:gen]],
        pdmax = isempty(ref[:storage]) ? Vector{NamedTuple{(:i,), Tuple{Int64}}}() : [s["charge_rating"] for (i, s) in ref[:storage]],
        pcmax = isempty(ref[:storage]) ? Vector{NamedTuple{(:i,), Tuple{Int64}}}() : [s["discharge_rating"] for (i, s) in ref[:storage]],
        srating = isempty(ref[:storage]) ? Vector{NamedTuple{(:i,), Tuple{Int64}}}() : [s["thermal_rating"] for (i, s) in ref[:storage]],
        emax = isempty(ref[:storage]) ? Vector{NamedTuple{(:i,), Tuple{Int64}}}() : [s["energy_rating"] for (i, s) in ref[:storage]],
    )

    @info "Saving JLD2 cache file"
    d, f = splitdir(filename)
    name,ext = splitext(f)
    JLD2.save(joinpath(TMPDIR, name * ".jld2"), "data", data, "dicts", dicts)

    return data, dicts
end

