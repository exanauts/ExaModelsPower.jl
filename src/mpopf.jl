using DelimitedFiles

function parse_mp_power_data(filename, N, corrective_action_ratio, pd, qd, backend)

    data = parse_ac_power_data(filename, nothing)
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
    
    update_load_data(data.busarray, pd, qd, data.baseMVA[])
    
    return convert_data(data,backend)
end

function update_load_data(busarray, pd, qd, baseMVA)
    for (i,b) in enumerate(busarray)
        busarray[i] = (;b..., pd = pd[b.bus_i, b.t] / baseMVA, qd = qd[b.bus_i, b.t] / baseMVA)
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
    
    w = ExaModels.ExaCore(T; backend = backend)

    va = ExaModels.variable(w, length(data.bus), N;)

    vm = ExaModels.variable(
        w,
        length(data.bus), N;
        start = repeat(fill!(similar(data.bus, Float64), 1.0), outer = (N,1)),
        lvar = repeat(data.vmin, outer = (N,1)),
        uvar = repeat(data.vmax, outer = (N,1)),
    )
    pg = ExaModels.variable(w, length(data.gen), N; lvar = repeat(data.pmin, outer = (N,1)), uvar = repeat(data.pmax, outer = (N,1)))

    qg = ExaModels.variable(w, length(data.gen), N; lvar = repeat(data.qmin, outer = (N,1)), uvar = repeat(data.qmax, outer = (N,1))) 

    p = ExaModels.variable(w, length(data.arc), N; lvar = repeat(-data.rate_a, outer = (N,1)), uvar = repeat(data.rate_a, outer = (N,1)))

    q = ExaModels.variable(w, length(data.arc), N; lvar = repeat(-data.rate_a, outer = (N,1)), uvar = repeat(data.rate_a, outer = (N,1)))

    

    o = ExaModels.objective(
        w,
        g.cost1 * pg[g.i,g.t]^2 + g.cost2 * pg[g.i,g.t] + g.cost3 for g in data.genarray
    )

    c1 = ExaModels.constraint(w, va[i,t] for (i,t) in data.refarray)

    c2 = ExaModels.constraint(
        w,
        p[b.f_idx, b.t] - b.c5 * vm[b.f_bus, b.t]^2 -
        b.c3 * (vm[b.f_bus, b.t] * vm[b.t_bus, b.t] * cos(va[b.f_bus, b.t] - va[b.t_bus, b.t])) -
        b.c4 * (vm[b.f_bus, b.t] * vm[b.t_bus, b.t] * sin(va[b.f_bus, b.t] - va[b.t_bus, b.t])) for
            b in data.barray
    )

    c3 = ExaModels.constraint(
        w,
        q[b.f_idx, b.t] +
        b.c6 * vm[b.f_bus, b.t]^2 +
        b.c4 * (vm[b.f_bus, b.t] * vm[b.t_bus, b.t] * cos(va[b.f_bus, b.t] - va[b.t_bus, b.t])) -
        b.c3 * (vm[b.f_bus, b.t] * vm[b.t_bus, b.t] * sin(va[b.f_bus, b.t] - va[b.t_bus, b.t])) for
        b in data.barray
    )

    c4 = ExaModels.constraint(
        w,
        p[b.t_idx, b.t] - b.c7 * vm[b.t_bus, b.t]^2 -
        b.c1 * (vm[b.t_bus, b.t] * vm[b.f_bus, b.t] * cos(va[b.t_bus, b.t] - va[b.f_bus, b.t])) -
        b.c2 * (vm[b.t_bus, b.t] * vm[b.f_bus, b.t] * sin(va[b.t_bus, b.t] - va[b.f_bus, b.t])) for
        b in data.barray
    )

    c5 = ExaModels.constraint(
        w,
        q[b.t_idx, b.t] +
        b.c8 * vm[b.t_bus, b.t]^2 +
        b.c2 * (vm[b.t_bus, b.t] * vm[b.f_bus, b.t] * cos(va[b.t_bus, b.t] - va[b.f_bus, b.t])) -
        b.c1 * (vm[b.t_bus, b.t] * vm[b.f_bus, b.t] * sin(va[b.t_bus, b.t] - va[b.f_bus, b.t])) for
        b in data.barray
    )

    c6 = ExaModels.constraint(
        w,
        va[b.f_bus, b.t] - va[b.t_bus, b.t] for b in data.barray;
        lcon = repeat(data.angmin, outer = (N,1)),
        ucon = repeat(data.angmax, outer = (N,1)),
    )
    c7 = ExaModels.constraint(
        w,
        p[b.f_idx, b.t]^2 + q[b.f_idx, b.t]^2 - b.rate_a_sq for b in data.barray;
        lcon = repeat(fill!(similar(data.branch, Float64, length(data.branch)), -Inf), outer = (N,1)),
    )
    c8 = ExaModels.constraint(
        w,
        p[b.t_idx, b.t]^2 + q[b.t_idx, b.t]^2 - b.rate_a_sq for b in data.barray;
        lcon = repeat(fill!(similar(data.branch, Float64, length(data.branch)), -Inf), outer = (N,1)),
    )

    c9 = ExaModels.constraint(w, b.pd + b.gs * vm[b.i, b.t]^2 for b in data.busarray)
    c10 = ExaModels.constraint(w, b.qd - b.bs * vm[b.i, b.t]^2 for b in data.busarray)

    c11 = ExaModels.constraint!(w, c9, a.cindex => p[a.i, a.t] for a in data.arcarray)
    c12 = ExaModels.constraint!(w, c10, a.cindex => q[a.i, a.t] for a in data.arcarray)

    c13 = ExaModels.constraint!(w, c9, g.cindex => -pg[g.i, g.t] for g in data.genarray)
    c14 = ExaModels.constraint!(w, c10, g.cindex => -qg[g.i, g.t] for g in data.genarray)

    
    c15 = ExaModels.constraint(
        w,
        pg[g.i, g.t1] - pg[g.i, g.t2] for g in data.idx;
        lcon = repeat(-data.Δp, outer = (N,1)),
        ucon = repeat( data.Δp, outer = (N,1)),
    )

    return ExaModels.ExaModel(w; kwargs...)

end
