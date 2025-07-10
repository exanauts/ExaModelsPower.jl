function build_polar_opf(data :: Data; backend = nothing, T=Float64, kwargs...)
    core = ExaCore(T; backend = backend)

    va = variable(core, length(data.bus);)
    vm = variable(
            core,
            length(data.bus);
            start = fill!(similar(data.bus, Float64), 1.0),
            lvar = [b.vmin for b in data.bus],
            uvar = [b.vmax for b in data.bus],
        )

    pg = variable(core, length(data.gen); lvar = [g.pmin for g in data.gen], uvar = [g.pmax for g in data.gen])
    qg = variable(core, length(data.gen); lvar = [g.qmin for g in data.gen], uvar = [g.qmax for g in data.gen])

    p = variable(core, length(data.arc); lvar = -data.ratea, uvar = data.ratea)
    q = variable(core, length(data.arc); lvar = -data.ratea, uvar = data.ratea)

    objective(core, gen_cost(g, pg[g.i]) for g in data.gen)

    c_ref_angle = constraint(core, c_ref_angle_polar(va[i]) for i in data.ref_buses)

    c_to_active_power_flow = constraint(core, c_to_active_power_flow_polar(b, p[b.fbus],
        vm[b.fbus], vm[b.tbus], va[b.fbus], va[b.tbus]) for b in data.branch)

    c_to_reactive_power_flow = constraint(core, c_to_reactive_power_flow_polar(b, q[b.fbus],
        vm[b.fbus],vm[b.tbus],va[b.fbus],va[b.tbus]) for b in data.branch)

    c_from_active_power_flow = constraint(core, c_from_active_power_flow_polar(b, p[b.tbus],
        vm[b.fbus],vm[b.tbus],va[b.fbus],va[b.tbus]) for b in data.branch)

    c_from_reactive_power_flow = constraint(core, c_from_reactive_power_flow_polar(b, q[b.tbus],
        vm[b.fbus],vm[b.tbus],va[b.fbus],va[b.tbus]) for b in data.branch)

    c_phase_angle_diff = constraint(
        core,
        c_phase_angle_diff_polar(b,va[b.f_bus],va[b.t_bus]) for b in data.branch;
        lcon = [b.angmin for b in data.branch],
        ucon = [b.angmax for b in data.branch],
    )

    c_active_power_balance = constraint(core, c_active_power_balance_demand_polar(b, vm[i]) for (i, b) in (data.bus))

    c_reactive_power_balance = constraint(core, c_reactive_power_balance_demand_polar(b, vm[i]) for (i, b) in enumerate(data.bus))

    # c_active_power_balance_arcs = constraint!(core, c_active_power_balance, a.bus => p[a.i] for a in data.arc)
    # c_reactive_power_balance_arcs = constraint!(core, c_reactive_power_balance, a.bus => q[a.i] for a in data.arc)

    c_active_power_balance_gen = constraint!(core, c_active_power_balance, g.bus => -pg[i] for (i, g) in enumerate(data.gen))
    c_active_power_balance_gen = constraint!(core, c_reactive_power_balance, g.bus => -qg[i] for (i, g) in enumerate(data.gen))

    c_from_thermal_limit = constraint(
        core, c_thermal_limit(b, p[b.fbus], q[b.fbus]) for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
        )

    c_to_thermal_limit = constraint(
        core, c_thermal_limit(b,p[b.tbus],q[b.tbus])
        for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    model = ExaModel(core; kwargs...)
    
    vars = (
            va = va,
            vm = vm,
            pg = pg,
            qg = qg,
            p = p,        
            q = q
        )

    cons = (
        c_ref_angle = c_ref_angle,
        c_to_active_power_flow = c_to_active_power_flow,
        c_to_reactive_power_flow = c_to_reactive_power_flow,
        c_from_active_power_flow = c_from_active_power_flow,
        c_from_reactive_power_flow = c_from_reactive_power_flow,
        c_phase_angle_diff = c_phase_angle_diff,
        c_active_power_balance = c_active_power_balance,
        c_reactive_power_balance = c_reactive_power_balance,
        c_from_thermal_limit = c_from_thermal_limit,
        c_to_thermal_limit = c_to_thermal_limit
    )

    return model, vars, cons
end

function build_rect_opf(data; backend = nothing, T=Float64, kwargs...)
    core = ExaCore(T; backend = backend)

    vr = variable(core, length(data.bus), start = fill!(similar(data.bus, Float64), 1.0))
    vim = variable(core, length(data.bus))

    pg = variable(core, length(data.gen); lvar = [g.pmin for g in data.gen], uvar = [g.pmax for g in data.gen])
    qg = variable(core, length(data.gen); lvar = [g.qmin for g in data.gen], uvar = [g.qmax for g in data.gen])

    p = variable(core, length(data.arc); lvar = -data.ratea, uvar = data.ratea)
    q = variable(core, length(data.arc); lvar = -data.ratea, uvar = data.ratea)

    o = objective(core, gen_cost(g, pg[i]) for (i, g) in enumerate(data.gen))

    c_ref_angle = constraint(core, c_ref_angle_rect(vr[i], vim[i]) for i in data.ref_buses)

    c_to_active_power_flow = constraint(core, c_to_active_power_flow_rect(b,p[b.f_idx],
        vr[b.f_bus],vr[b.t_bus],vim[b.f_bus],vim[b.t_bus]) for b in data.branch)

    c_to_reactive_power_flow = constraint(core, c_to_reactive_power_flow_rect(b,q[b.f_idx],
        vr[b.f_bus],vr[b.t_bus],vim[b.f_bus],vim[b.t_bus]) for b in data.branch)

    c_from_active_power_flow = constraint(core, c_from_active_power_flow_rect(b,p[b.t_idx],
        vr[b.f_bus],vr[b.t_bus],vim[b.f_bus],vim[b.t_bus]) for b in data.branch)
    
    c_from_reactive_power_flow = constraint(core, c_from_reactive_power_flow_rect(b,q[b.t_idx],
        vr[b.f_bus],vr[b.t_bus],vim[b.f_bus],vim[b.t_bus]) for b in data.branch)

    c_phase_angle_diff = constraint(
        core, c_phase_angle_diff_rect(b,
        vr[b.f_bus],vr[b.t_bus],vim[b.f_bus],vim[b.t_bus])
        for b in data.branch;
        lcon = [b.angmin for b in data.branch],
        ucon = [b.angmax for b in data.branch],
    )

    c_active_power_balance = constraint(core, c_active_power_balance_demand_rect(b, vr[b.i], vim[b.i]) for b in data.bus)

    c_reactive_power_balance = constraint(core, c_reactive_power_balance_demand_rect(b, vr[b.i], vim[b.i]) for b in data.bus)

    # c_active_power_balance_arcs = constraint!(core, c_active_power_balance, a.bus => p[a.i] for a in data.arc)
    # c_reactive_power_balance_arcs = constraint!(core, c_reactive_power_balance, a.bus => q[a.i] for a in data.arc)

    # c_active_power_balance_gen = constraint!(core, c_active_power_balance, g.bus => -pg[i] for (i, g) in enumerate(data.gen))
    # c_reactive_power_balance_gen = constraint!(core, c_reactive_power_balance, g.bus => -qg[i] for (i, g) in enumerate(data.gen))

    c_from_thermal_limit = constraint(
        core, c_thermal_limit(b,p[b.fbus], q[b.fbus]) for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c_to_thermal_limit = constraint(
        core, c_thermal_limit(b,p[b.tbus], q[b.tbus])
        for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c_voltage_magnitude = constraint(
        core, c_voltage_magnitude_rect(vr[b.i], vim[b.i]) for b in data.bus; 
        lcon = [b.vmin for b in data.bus].^2,
        ucon = [b.vmax for b in data.bus].^2,
    )
    
    model = ExaModel(core; kwargs...)
    
    vars = (
        vr = vr,
        vim = vim,
        pg = pg,
        qg = qg,
        p = p,        
        q = q
    )

    cons = (
        c_ref_angle = c_ref_angle,
        c_to_active_power_flow = c_to_active_power_flow,
        c_to_reactive_power_flow = c_to_reactive_power_flow,
        c_from_active_power_flow = c_from_active_power_flow,
        c_from_reactive_power_flow = c_from_reactive_power_flow,
        c_phase_angle_diff = c_phase_angle_diff,
        c_active_power_balance = c_active_power_balance,
        c_reactive_power_balance = c_reactive_power_balance,
        c_from_thermal_limit = c_from_thermal_limit,
        c_to_thermal_limit = c_to_thermal_limit,
        c_voltage_magnitude = c_voltage_magnitude
    )

    return model, vars, cons
end

function opf_model(
    filename;
    backend = nothing,
    T = Float64,
    form = :polar,
    kwargs...,
)
    data :: Data{T} = parse_ac_power_data(T, filename)

    if form == :polar
        return build_polar_opf(data, backend = backend, T=T, kwargs...)
    elseif form == :rect
        return build_rect_opf(data, backend = backend, T=T, kwargs...)
    else
        error("Invalid coordinate symbol - valid options are :polar or :rect")
    end
end
