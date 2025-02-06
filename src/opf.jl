function opf_model(
    filename;
    backend = nothing,
    T = Float64,
    coords = "polar",
    kwargs...,
)

    data, _ = parse_ac_power_data(filename)
    data = convert_data(data, backend)
    
    core = ExaModels.ExaCore(T; backend = backend)

    if coords == "polar"
        va = ExaModels.variable(core, length(data.bus);)

        vm = ExaModels.variable(
            core,
            length(data.bus);
            start = fill!(similar(data.bus, Float64), 1.0),
            lvar = data.vmin,
            uvar = data.vmax,
        )
    elseif coords == "rect"
        vr = ExaModels.variable(core, length(data.bus), start = fill!(similar(data.bus, Float64), 1.0))

        vim = ExaModels.variable(core, length(data.bus))
    else
        error("Invalid coordinate system - valid options are 'polar' or 'rect'")
    end

    pg = ExaModels.variable(core, length(data.gen); lvar = data.pmin, uvar = data.pmax)

    qg = ExaModels.variable(core, length(data.gen); lvar = data.qmin, uvar = data.qmax)

    p = ExaModels.variable(core, length(data.arc); lvar = -data.rate_a, uvar = data.rate_a)

    q = ExaModels.variable(core, length(data.arc); lvar = -data.rate_a, uvar = data.rate_a)

    o = ExaModels.objective(
        core,
        g.cost1 * pg[g.i]^2 + g.cost2 * pg[g.i] + g.cost3 for g in data.gen
    )

    if coords == "polar"
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

        c7 = ExaModels.constraint(core, b.pd + b.gs * vm[b.i]^2 for b in data.bus)

        c8 = ExaModels.constraint(core, b.qd - b.bs * vm[b.i]^2 for b in data.bus)

    elseif coords == "rect"
        c1 = ExaModels.constraint(core, atan(vim[i]/vr[i]) for i in data.ref_buses)

        c2 = ExaModels.constraint(
            core,
            p[b.f_idx] - b.c5 * (vr[b.f_bus]^2+vim[b.f_bus]^2) -
            b.c3 * (vr[b.f_bus]*vr[b.t_bus] + vim[b.f_bus]*vim[b.t_bus]) -
            b.c4 * (vim[b.f_bus]*vr[b.t_bus] - vr[b.f_bus]*vim[b.t_bus]) for
            b in data.branch
        )

        c3 = ExaModels.constraint(
            core,
            q[b.f_idx] +
            b.c6 * (vr[b.f_bus]^2+vim[b.f_bus]^2) +
            b.c4 * (vr[b.f_bus]*vr[b.t_bus] + vim[b.f_bus]*vim[b.t_bus]) -
            b.c3 * (vim[b.f_bus]*vr[b.t_bus] - vr[b.f_bus]*vim[b.t_bus]) for
            b in data.branch
        )

        c4 = ExaModels.constraint(
            core,
            p[b.t_idx] - b.c7 * (vr[b.t_bus]^2+vim[b.t_bus]^2) -
            b.c1 * (vr[b.f_bus]*vr[b.t_bus] + vim[b.f_bus]*vim[b.t_bus]) -
            b.c2 * (vim[b.t_bus]*vr[b.f_bus] - vr[b.t_bus]*vim[b.f_bus]) for
            b in data.branch
        )

        c5 = ExaModels.constraint(
            core,
            q[b.t_idx] +
            b.c8 * (vr[b.t_bus]^2+vim[b.t_bus]^2) +
            b.c2 * (vr[b.f_bus]*vr[b.t_bus] + vim[b.f_bus]*vim[b.t_bus]) -
            b.c1 * (vim[b.t_bus]*vr[b.f_bus] - vr[b.t_bus]*vim[b.f_bus]) for
            b in data.branch
        )

        c6 = ExaModels.constraint(
            core,
            atan(vim[b.f_bus]/vr[b.f_bus]) - atan(vim[b.t_bus]/vr[b.t_bus]) for b in data.branch;
            lcon = data.angmin,
            ucon = data.angmax,
        )

        c7 = ExaModels.constraint(core, b.pd + b.gs * (vr[b.i]^2 + vim[b.i]^2) for b in data.bus)

        c8 = ExaModels.constraint(core, b.qd - b.bs * (vr[b.i]^2 + vim[b.i]^2) for b in data.bus)

        c11 = ExaModels.constraint(
            core, 
            vr[b.i]^2 + vim[b.i]^2 for b in data.bus; 
            lcon = data.vmin.^2, 
            ucon = data.vmax.^2
        ) 
    end

    c7a = ExaModels.constraint!(core, c7, a.bus => p[a.i] for a in data.arc)
    c8a = ExaModels.constraint!(core, c8, a.bus => q[a.i] for a in data.arc)

    c7b = ExaModels.constraint!(core, c7, g.bus => -pg[g.i] for g in data.gen)
    c8b = ExaModels.constraint!(core, c8, g.bus => -qg[g.i] for g in data.gen)

    c9 = ExaModels.constraint(
        core,
        p[b.f_idx]^2 + q[b.f_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )
    c10 = ExaModels.constraint(
        core,
        p[b.t_idx]^2 + q[b.t_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    model =ExaModels.ExaModel(core; kwargs...)
    
    if coords == "polar"
        vars = (
            va = va,
            vm = vm,
            pg = pg,
            qg = qg,
            p = p,        
            q = q
        )
    elseif coords == "rect"
        vars = (
            vr = vr,
            vim = vim,
            pg = pg,
            qg = qg,
            p = p,        
            q = q
        )
    end

    return model, vars
end
