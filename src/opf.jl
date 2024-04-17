function opf_model(
    filename;
    backend = nothing,
    T = Float64,
    kwargs...,
)

    data, _ = parse_ac_power_data(filename)
    data = convert_data(data, backend)
    
    core = ExaModels.ExaCore(T; backend = backend)

    va = ExaModels.variable(core, length(data.bus);)

    vm = ExaModels.variable(
        core,
        length(data.bus);
        start = fill!(similar(data.bus, Float64), 1.0),
        lvar = data.vmin,
        uvar = data.vmax,
    )
    pg = ExaModels.variable(core, length(data.gen); lvar = data.pmin, uvar = data.pmax)

    qg = ExaModels.variable(core, length(data.gen); lvar = data.qmin, uvar = data.qmax)

    p = ExaModels.variable(core, length(data.arc); lvar = -data.rate_a, uvar = data.rate_a)

    q = ExaModels.variable(core, length(data.arc); lvar = -data.rate_a, uvar = data.rate_a)

    o = ExaModels.objective(
        core,
        g.cost1 * pg[g.i]^2 + g.cost2 * pg[g.i] + g.cost3 for g in data.gen
    )

    c1 = ExaModels.constraint(core, va[i] for i in data.ref_buses)

    c2 = ExaModels.constraint(
        core,
        p[b.f_idx] - b.c5 * vm[b.f_bus]^2 -
        b.c3 * (vm[b.f_bus] * vm[b.t_bus] * cos(va[b.f_bus] - va[b.t_bus])) -
        b.c4 * (vm[b.f_bus] * vm[b.t_bus] * sin(va[b.f_bus] - va[b.t_bus])) for
        b in data.branch
    )

    c3 = ExaModels.constraint(
        core,
        q[b.f_idx] +
        b.c6 * vm[b.f_bus]^2 +
        b.c4 * (vm[b.f_bus] * vm[b.t_bus] * cos(va[b.f_bus] - va[b.t_bus])) -
        b.c3 * (vm[b.f_bus] * vm[b.t_bus] * sin(va[b.f_bus] - va[b.t_bus])) for
        b in data.branch
    )

    c4 = ExaModels.constraint(
        core,
        p[b.t_idx] - b.c7 * vm[b.t_bus]^2 -
        b.c1 * (vm[b.t_bus] * vm[b.f_bus] * cos(va[b.t_bus] - va[b.f_bus])) -
        b.c2 * (vm[b.t_bus] * vm[b.f_bus] * sin(va[b.t_bus] - va[b.f_bus])) for
        b in data.branch
    )

    c5 = ExaModels.constraint(
        core,
        q[b.t_idx] +
        b.c8 * vm[b.t_bus]^2 +
        b.c2 * (vm[b.t_bus] * vm[b.f_bus] * cos(va[b.t_bus] - va[b.f_bus])) -
        b.c1 * (vm[b.t_bus] * vm[b.f_bus] * sin(va[b.t_bus] - va[b.f_bus])) for
        b in data.branch
    )

    c6 = ExaModels.constraint(
        core,
        va[b.f_bus] - va[b.t_bus] for b in data.branch;
        lcon = data.angmin,
        ucon = data.angmax,
    )
    c7 = ExaModels.constraint(
        core,
        p[b.f_idx]^2 + q[b.f_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )
    c8 = ExaModels.constraint(
        core,
        p[b.t_idx]^2 + q[b.t_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c9 = ExaModels.constraint(core, b.pd + b.gs * vm[b.i]^2 for b in data.bus)

    c10 = ExaModels.constraint(core, b.qd - b.bs * vm[b.i]^2 for b in data.bus)

    c11 = ExaModels.constraint!(core, c9, a.bus => p[a.i] for a in data.arc)
    c12 = ExaModels.constraint!(core, c10, a.bus => q[a.i] for a in data.arc)

    c13 = ExaModels.constraint!(core, c9, g.bus => -pg[g.i] for g in data.gen)
    c14 = ExaModels.constraint!(core, c10, g.bus => -qg[g.i] for g in data.gen)

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
