using DelimitedFiles

function parse_mp_power_data(filename, N, corrective_action_ratio, bus, pd, qd, backend)

    data, dicts = parse_ac_power_data(filename)

    nbus = length(data.bus)

    @assert nbus == size(pd,1)
    
    data = (
        ;
        data...,
        refarray = [(i,t) for t in 1:N for i in data.ref_buses],
        barray = [(;b..., t = t) for t in 1:N for b in data.branch],
        busarray = [(;b..., t = t) for t in 1:N for b in data.bus],
        arcarray = [(;a..., t = t, cindex = (t-1) * nbus + a.bus) for t in 1:N for a in data.arc],
        genarray = [(;g..., t = t, cindex = (t-1) * nbus + g.bus) for t in 1:N for g in data.gen],
        idx = [(t1=t, t2=t+1, i = g.i) for t in 1:N-1 for g in data.gen],
        Δp = corrective_action_ratio .* (data.pmax .- data.pmin)
    )

    busmap = zeros(Int, nbus)

    for j=1:length(busmap)
        busmap[dicts.bus[Int(bus[j,1])]] = j
    end
    
    update_load_data(data.busarray, pd, qd, data.baseMVA[], busmap)
    
    return convert_data(data,backend)
end

function update_load_data(busarray, pd, qd, baseMVA, busmap)
    for (i,b) in enumerate(busarray)
        busarray[i] = (;b..., pd = pd[busmap[b.i], b.t] / baseMVA, qd = qd[busmap[b.i], b.t] / baseMVA)
    end
end


function mpopf_model(
    filename, bus_data, active_power_data, reactive_power_data;
    bus = readdlm(bus_data),
    pd = readdlm(active_power_data),
    qd = readdlm(reactive_power_data),
    N = size(pd,2), 
    corrective_action_ratio = 0.1,
    backend = nothing,
    T = Float64,
    kwargs...,
)


    data = parse_mp_power_data(filename, N, corrective_action_ratio, bus, pd, qd, backend)
    
    core = ExaModels.ExaCore(T; backend = backend)

    va = ExaModels.variable(core, length(data.bus), N;)

    vm = ExaModels.variable(
        core,
        length(data.bus), N;
        start = repeat(fill!(similar(data.bus, Float64), 1.0), outer = (N,1)),
        lvar = repeat(data.vmin, outer = (N,1)),
        uvar = repeat(data.vmax, outer = (N,1)),
    )
    pg = ExaModels.variable(core, length(data.gen), N; lvar = repeat(data.pmin, outer = (N,1)), uvar = repeat(data.pmax, outer = (N,1)))

    qg = ExaModels.variable(core, length(data.gen), N; lvar = repeat(data.qmin, outer = (N,1)), uvar = repeat(data.qmax, outer = (N,1))) 

    p = ExaModels.variable(core, length(data.arc), N; lvar = repeat(-data.rate_a, outer = (N,1)), uvar = repeat(data.rate_a, outer = (N,1)))

    q = ExaModels.variable(core, length(data.arc), N; lvar = repeat(-data.rate_a, outer = (N,1)), uvar = repeat(data.rate_a, outer = (N,1)))

    

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
        lcon = repeat(data.angmin, outer = (N,1)),
        ucon = repeat(data.angmax, outer = (N,1)),
    )
    c7 = ExaModels.constraint(
        core,
        p[b.f_idx, b.t]^2 + q[b.f_idx, b.t]^2 - b.rate_a_sq for b in data.barray;
        lcon = repeat(fill!(similar(data.branch, Float64, length(data.branch)), -Inf), outer = (N,1)),
    )
    c8 = ExaModels.constraint(
        core,
        p[b.t_idx, b.t]^2 + q[b.t_idx, b.t]^2 - b.rate_a_sq for b in data.barray;
        lcon = repeat(fill!(similar(data.branch, Float64, length(data.branch)), -Inf), outer = (N,1)),
    )

    c9 = ExaModels.constraint(core, b.pd + b.gs * vm[b.i, b.t]^2 for b in data.busarray)
    c10 = ExaModels.constraint(core, b.qd - b.bs * vm[b.i, b.t]^2 for b in data.busarray)

    c11 = ExaModels.constraint!(core, c9, a.cindex => p[a.i, a.t] for a in data.arcarray)
    c12 = ExaModels.constraint!(core, c10, a.cindex => q[a.i, a.t] for a in data.arcarray)

    c13 = ExaModels.constraint!(core, c9, g.cindex => -pg[g.i, g.t] for g in data.genarray)
    c14 = ExaModels.constraint!(core, c10, g.cindex => -qg[g.i, g.t] for g in data.genarray)

    
    c15 = ExaModels.constraint(
        core,
        pg[g.i, g.t1] - pg[g.i, g.t2] for g in data.idx;
        lcon = repeat(-data.Δp, outer = (N,1)),
        ucon = repeat( data.Δp, outer = (N,1)),
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
