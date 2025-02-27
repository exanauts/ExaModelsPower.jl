function build_polar_opf(data; backend = nothing, T=Float64, kwargs...)
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
        core, obj(g, pg[g.i]) for g in data.gen)
    c1 = ExaModels.constraint(core, c1_polar(va[i]) for i in data.ref_buses)

    c2 = ExaModels.constraint(core, c2_polar(b, p[b.f_idx],
        vm[b.f_bus],vm[b.t_bus],va[b.f_bus],va[b.t_bus]) for b in data.branch)

    c3 = ExaModels.constraint(core, c3_polar(b, q[b.f_idx],
        vm[b.f_bus],vm[b.t_bus],va[b.f_bus],va[b.t_bus]) for b in data.branch)

    c4 = ExaModels.constraint(core, c4_polar(b, p[b.t_idx],
        vm[b.f_bus],vm[b.t_bus],va[b.f_bus],va[b.t_bus]) for b in data.branch)

    c5 = ExaModels.constraint(core, c5_polar(b, q[b.t_idx],
        vm[b.f_bus],vm[b.t_bus],va[b.f_bus],va[b.t_bus]) for b in data.branch)

    c6 = ExaModels.constraint(
        core,
        c6_polar(b,va[b.f_bus],va[b.t_bus]) for b in data.branch;
        lcon = data.angmin,
        ucon = data.angmax,
    )

    c7 = ExaModels.constraint(core, c7_polar(b, vm[b.i]) for b in data.bus)

    c8 = ExaModels.constraint(core, c8_polar(b, vm[b.i]) for b in data.bus)

    c7a = ExaModels.constraint!(core, c7, a.bus => p[a.i] for a in data.arc)
    c8a = ExaModels.constraint!(core, c8, a.bus => q[a.i] for a in data.arc)

    c7b = ExaModels.constraint!(core, c7, g.bus => -pg[g.i] for g in data.gen)
    c8b = ExaModels.constraint!(core, c8, g.bus => -qg[g.i] for g in data.gen)

    c9 = ExaModels.constraint(
        core, c9_10(b,p[b.f_idx],q[b.f_idx]) for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
        )

    c10 = ExaModels.constraint(
        core, c9_10(b,p[b.t_idx],q[b.t_idx])
        for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
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

function build_rect_opf(data; backend = nothing, T=Float64, kwargs...)
    core = ExaModels.ExaCore(T; backend = backend)

    vr = ExaModels.variable(core, length(data.bus), start = fill!(similar(data.bus, Float64), 1.0))

    vim = ExaModels.variable(core, length(data.bus))

    pg = ExaModels.variable(core, length(data.gen); lvar = data.pmin, uvar = data.pmax)

    qg = ExaModels.variable(core, length(data.gen); lvar = data.qmin, uvar = data.qmax)

    p = ExaModels.variable(core, length(data.arc); lvar = -data.rate_a, uvar = data.rate_a)

    q = ExaModels.variable(core, length(data.arc); lvar = -data.rate_a, uvar = data.rate_a)

    o = ExaModels.objective(
        core, obj(g, pg[g.i]) for g in data.gen)

    c1 = ExaModels.constraint(core, c1_rect(vr[i], vim[i]) for i in data.ref_buses)

    c2 = ExaModels.constraint(core, c2_rect(b,p[b.f_idx],
        vr[b.f_bus],vr[b.t_bus],vim[b.f_bus],vim[b.t_bus]) for b in data.branch)

    c3 = ExaModels.constraint(core, c3_rect(b,q[b.f_idx],
        vr[b.f_bus],vr[b.t_bus],vim[b.f_bus],vim[b.t_bus]) for b in data.branch)

    c4 = ExaModels.constraint(core, c4_rect(b,p[b.t_idx],
        vr[b.f_bus],vr[b.t_bus],vim[b.f_bus],vim[b.t_bus]) for b in data.branch)
    
    c5 = ExaModels.constraint(core, c5_rect(b,q[b.t_idx],
        vr[b.f_bus],vr[b.t_bus],vim[b.f_bus],vim[b.t_bus]) for b in data.branch)

    c6 = ExaModels.constraint(
        core, c6_rect(b,
        vr[b.f_bus],vr[b.t_bus],vim[b.f_bus],vim[b.t_bus])
        for b in data.branch;
        lcon = data.angmin,
        ucon = data.angmax,
    )

    c7 = ExaModels.constraint(core, c7_rect(b, vr[b.i], vim[b.i]) for b in data.bus)

    c8 = ExaModels.constraint(core, c8_rect(b, vr[b.i], vim[b.i]) for b in data.bus)

    c7a = ExaModels.constraint!(core, c7, a.bus => p[a.i] for a in data.arc)
    c8a = ExaModels.constraint!(core, c8, a.bus => q[a.i] for a in data.arc)

    c7b = ExaModels.constraint!(core, c7, g.bus => -pg[g.i] for g in data.gen)
    c8b = ExaModels.constraint!(core, c8, g.bus => -qg[g.i] for g in data.gen)

    c9 = ExaModels.constraint(
        core, c9_10(b,p[b.f_idx], q[b.f_idx]) for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c10 = ExaModels.constraint(
        core, c9_10(b,p[b.t_idx], q[b.t_idx])
        for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c11 = ExaModels.constraint(
        core, c11_rect(vr[b.i], vim[b.i]) for b in data.bus; 
        lcon = data.vmin.^2, 
        ucon = data.vmax.^2
    ) 
    
    model =ExaModels.ExaModel(core; kwargs...)
    
    vars = (
        vr = vr,
        vim = vim,
        pg = pg,
        qg = qg,
        p = p,        
        q = q
    )

    return model, vars
end

function opf_model(
    filename;
    backend = nothing,
    T = Float64,
    symbol = "polar",
    kwargs...,
)

    data, _ = parse_ac_power_data(filename)
    data = convert_data(data, backend)

    if symbol == "polar"
        return build_polar_opf(data, backend = backend, T=T, kwargs...)
    elseif symbol == "rect"
        return build_rect_opf(data, backend = backend, T=T, kwargs...)
    else
        error("Invalid coordinate symbol - valid options are 'polar' or 'rect'")
    end
end