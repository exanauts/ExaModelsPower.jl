using DelimitedFiles

function parse_mp_power_data(filename, N, corrective_action_ratio, pd, qd, backend)

    data, dicts = parse_ac_power_data(filename)

    nbus = length(data.bus)

    @assert nbus == size(pd,1)
    
    data = (
        ;
        data...,
        refarray = [(i,t) for i in data.ref_buses, t in 1:N],
        barray = [(;b..., t = t) for b in data.branch, t in 1:N ],
        busarray = [(;b..., t = t) for b in data.bus, t in 1:N ],
        arcarray = [(;a..., t = t) for a in data.arc, t in 1:N ],
        genarray = [(;g..., t = t) for g in data.gen, t in 1:N ],
        Δp = corrective_action_ratio .* (data.pmax .- data.pmin)
    )
    
    update_load_data(data.busarray, pd, qd, data.baseMVA[], dicts.bus)
    
    return convert_data(data,backend)
end

function update_load_data(busarray, pd, qd, baseMVA, busdict)
    for (idx ,pd_t) in pairs(pd)
        b = busarray[busdict[idx[1]], idx[2]]
        busarray[busdict[idx[1]], idx[2]] = (
                i = b.i,
                pd = pd_t/ baseMVA, 
                gs = b.gs,
                qd = qd[idx[1], idx[2]] / baseMVA,
                bs = b.bs,
                bus_type = b.bus_type,
                t = idx[2]
                )
    end
end


function mpopf_model(
    filename, active_power_data, reactive_power_data;
    pd = readdlm(active_power_data),
    qd = readdlm(reactive_power_data),
    N = size(pd,2), 
    corrective_action_ratio = 0.1,
    backend = nothing,
    T = Float64,
    kwargs...,
)


    data = parse_mp_power_data(filename, N, corrective_action_ratio, pd, qd, backend)
    
    core = ExaModels.ExaCore(T; backend = backend)

    va = ExaModels.variable(core, size(data.bus, 1), N;)

    vm = ExaModels.variable(
        core,
        size(data.bus, 1), N;
        start = ones(size(data.bus)),
        lvar = repeat(data.vmin, 1, N),
        uvar = repeat(data.vmax, 1, N),
    )
    pg = ExaModels.variable(core, size(data.gen, 1), N; lvar = repeat(data.pmin, 1, N), uvar = repeat(data.pmax, 1, N))

    qg = ExaModels.variable(core, size(data.gen, 1), N; lvar = repeat(data.qmin, 1, N), uvar = repeat(data.qmax, 1, N)) 

    p = ExaModels.variable(core, size(data.arc, 1), N; lvar = repeat(-data.rate_a, 1, N), uvar = repeat(data.rate_a, 1, N))

    q = ExaModels.variable(core, size(data.arc, 1), N; lvar = repeat(-data.rate_a, 1, N), uvar = repeat(data.rate_a, 1, N))

    

    o = ExaModels.objective(
        core,
        g.cost1 * pg[g.i,g.t]^2 + g.cost2 * pg[g.i,g.t] + g.cost3 for g in data.genarray
    )

    c1 = ExaModels.constraint(core, va[i,t] for (i,t) in data.refarray)

    c2 = ExaModels.constraint(
        core,
        p[b.f_idx, b.t] - b.c5 * vm[b.f_bus, b.t]^2 -
        b.c3 * (vm[b.f_bus, b.t] * vm[b.t_bus, b.t] * cos(va[b.f_bus, b.t] - va[b.t_bus, b.t])) -
        b.c4 * (vm[b.f_bus, b.t] * vm[b.t_bus, b.t] * sin(va[b.f_bus, b.t] - va[b.t_bus, b.t])) for
            b in data.barray
    )

    c3 = ExaModels.constraint(
        core,
        q[b.f_idx, b.t] +
        b.c6 * vm[b.f_bus, b.t]^2 +
        b.c4 * (vm[b.f_bus, b.t] * vm[b.t_bus, b.t] * cos(va[b.f_bus, b.t] - va[b.t_bus, b.t])) -
        b.c3 * (vm[b.f_bus, b.t] * vm[b.t_bus, b.t] * sin(va[b.f_bus, b.t] - va[b.t_bus, b.t])) for
        b in data.barray
    )

    c4 = ExaModels.constraint(
        core,
        p[b.t_idx, b.t] - b.c7 * vm[b.t_bus, b.t]^2 -
        b.c1 * (vm[b.t_bus, b.t] * vm[b.f_bus, b.t] * cos(va[b.t_bus, b.t] - va[b.f_bus, b.t])) -
        b.c2 * (vm[b.t_bus, b.t] * vm[b.f_bus, b.t] * sin(va[b.t_bus, b.t] - va[b.f_bus, b.t])) for
        b in data.barray
    )

    c5 = ExaModels.constraint(
        core,
        q[b.t_idx, b.t] +
        b.c8 * vm[b.t_bus, b.t]^2 +
        b.c2 * (vm[b.t_bus, b.t] * vm[b.f_bus, b.t] * cos(va[b.t_bus, b.t] - va[b.f_bus, b.t])) -
        b.c1 * (vm[b.t_bus, b.t] * vm[b.f_bus, b.t] * sin(va[b.t_bus, b.t] - va[b.f_bus, b.t])) for
        b in data.barray
    )

    c6 = ExaModels.constraint(
        core,
        va[b.f_bus, b.t] - va[b.t_bus, b.t] for b in data.barray;
        lcon = repeat(data.angmin, 1, N),
        ucon = repeat(data.angmax, 1, N),
    )
    

    c7 = ExaModels.constraint(core, b.pd + b.gs * vm[b.i, b.t]^2 for b in data.busarray)
    c8 = ExaModels.constraint(core, b.qd - b.bs * vm[b.i, b.t]^2 for b in data.busarray)

    c7a = ExaModels.constraint!(core, c7, a.bus + N*(a.t-1) => p[a.i, a.t] for a in data.arcarray)
    c8a = ExaModels.constraint!(core, c8, a.bus + N*(a.t-1) => q[a.i, a.t] for a in data.arcarray)

    c7b = ExaModels.constraint!(core, c7, g.bus + N*(g.t-1) => -pg[g.i, g.t] for g in data.genarray)
    c8b = ExaModels.constraint!(core, c8, g.bus + N*(g.t-1) => -qg[g.i, g.t] for g in data.genarray)

    c9 = ExaModels.constraint(
        core,
        p[b.f_idx, b.t]^2 + q[b.f_idx, b.t]^2 - b.rate_a_sq for b in data.barray;
        lcon = fill(-Inf, size(data.barray))
    )
    c10 = ExaModels.constraint(
        core,
        p[b.t_idx, b.t]^2 + q[b.t_idx, b.t]^2 - b.rate_a_sq for b in data.barray;
        lcon = fill(-Inf, size(data.barray))
    )

    c11 = ExaModels.constraint(
        core,
        pg[g.i, g.t -1] - pg[g.i, g.t] for g in data.genarray[:, 2:N];
        lcon = repeat(-data.Δp,  1, N-1),
        ucon = repeat( data.Δp, 1, N-1),
    )

    model =ExaModels.ExaModel(core; kwargs...)
    
    vars = (
        va = va,
        vm = vm,
        pg = pg,
        qg = qg,
        p = p,        
        q = q
    )

    return model, vars
end
