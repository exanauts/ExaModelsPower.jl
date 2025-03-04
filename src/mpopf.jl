using DelimitedFiles

#Curve as input
function parse_mp_power_data(filename, N, corrective_action_ratio, backend, curve)

    data, dicts = parse_ac_power_data(filename)

    nbus = length(data.bus)

    @assert length(curve) > 0

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
    
    update_load_data(data.busarray, curve)
    return convert_data(data,backend)
end

#Pd, Qd as inputs
function parse_mp_power_data(filename, N, corrective_action_ratio, pd, qd, backend)

    data, dicts = parse_ac_power_data(filename)

    nbus = length(data.bus)

    @assert nbus == size(pd, 1)

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

    return convert_data(data, backend)
end

#Curve as input
function update_load_data(busarray, curve)

    for t in eachindex(curve)
        for x in 1:size(busarray, 1)
            b = busarray[x, t]
            busarray[x, t] = (
                i = b.i,
                pd = b.pd*curve[t], 
                gs = b.gs,
                qd = b.qd*curve[t],
                bs = b.bs,
                bus_type = b.bus_type,
                t = t
                )
        end
    end
end

#Pd, Qd as input
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

function build_polar_mpopf(data, Nbus, N; backend = nothing, T = Float64, kwargs...)
    core = ExaCore(T; backend = backend)

    va = variable(core, Nbus, N;)
    vm = variable(
        core,
        Nbus, N;
        start = ones(size(data.busarray)),
        lvar = repeat(data.vmin, 1, N),
        uvar = repeat(data.vmax, 1, N),
    )

    pg = variable(core, size(data.gen, 1), N; lvar = repeat(data.pmin, 1, N), uvar = repeat(data.pmax, 1, N))
    qg = variable(core, size(data.gen, 1), N; lvar = repeat(data.qmin, 1, N), uvar = repeat(data.qmax, 1, N)) 

    p = variable(core, size(data.arc, 1), N; lvar = repeat(-data.rate_a, 1, N), uvar = repeat(data.rate_a, 1, N))
    q = variable(core, size(data.arc, 1), N; lvar = repeat(-data.rate_a, 1, N), uvar = repeat(data.rate_a, 1, N))

    o = objective(core, gen_cost(g, pg[g.i]) for g in data.genarray)

    c_ref_angle = constraint(core, c_ref_angle_polar(va[i,t]) for (i,t) in data.refarray)

    c_to_active_power_flow = constraint(core, c_to_active_power_flow_polar(b, p[b.f_idx, b.t],
        vm[b.f_bus, b.t], vm[b.t_bus, b.t], va[b.f_bus, b.t], va[b.t_bus, b.t]) for b in data.barray)
    
    c_to_reactive_power_flow = constraint(core, c_to_reactive_power_flow_polar(b, q[b.f_idx, b.t],
        vm[b.f_bus, b.t], vm[b.t_bus, b.t], va[b.f_bus, b.t], va[b.t_bus, b.t]) for b in data.barray)

    c_from_active_power_flow = constraint(core, c_from_active_power_flow_polar(b, p[b.t_idx, b.t],
        vm[b.f_bus, b.t], vm[b.t_bus, b.t], va[b.f_bus, b.t], va[b.t_bus, b.t]) for b in data.barray)
    
    c_from_reactive_power_flow = constraint(core, c_from_reactive_power_flow_polar(b, q[b.t_idx, b.t],
        vm[b.f_bus, b.t], vm[b.t_bus, b.t], va[b.f_bus, b.t], va[b.t_bus, b.t]) for b in data.barray)
    
    c_phase_angle_diff = constraint(
        core,
        c_phase_angle_diff_polar(b, va[b.f_bus, b.t], va[b.t_bus, b.t]) for b in data.barray;
        lcon = repeat(data.angmin, 1, N),
        ucon = repeat(data.angmax, 1, N),
    )
    
    c_active_power_balance = constraint(core, c_active_power_balance_demand_polar(b, vm[b.i, b.t]) for b in data.busarray)

    c_reactive_power_balance = constraint(core, c_reactive_power_balance_demand_polar(b, vm[b.i, b.t]) for b in data.busarray)

    c_active_power_balance_arcs = constraint!(core, c_active_power_balance, a.bus + Nbus*(a.t-1) => p[a.i, a.t] for a in data.arcarray)
    c_reactive_power_balance_arcs = constraint!(core, c_reactive_power_balance, a.bus + Nbus*(a.t-1) => q[a.i, a.t] for a in data.arcarray)

    c_active_power_balance_gen = constraint!(core, c_active_power_balance, g.bus + Nbus*(g.t-1) => -pg[g.i, g.t] for g in data.genarray)
    c_reactive_power_balance = constraint!(core, c_reactive_power_balance, g.bus + Nbus*(g.t-1) => -qg[g.i, g.t] for g in data.genarray)
    
    c_from_thermal_limit = constraint(
        core,
        c_thermal_limit(b, p[b.f_idx, b.t], q[b.f_idx, b.t]) for b in data.barray;
        lcon = fill(-Inf, size(data.barray))
    )

    c_to_thermal_limit = constraint(
        core,
        c_thermal_limit(b, p[b.t_idx, b.t], q[b.t_idx, b.t]) for b in data.barray;
        lcon = fill(-Inf, size(data.barray))
    )

    c_ramp_rate = constraint(
        core,
        c_ramp(pg[g.i, g.t -1], pg[g.i, g.t]) for g in data.genarray[:, 2:N];
        lcon = repeat(-data.Δp,  1, N-1),
        ucon = repeat( data.Δp, 1, N-1),
    )

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
            c_to_thermal_limit = c_to_thermal_limit,
            c_ramp_rate = c_ramp_rate
        )

    model =ExaModel(core; kwargs...)
    return model, vars, cons
end
    
function build_rect_mpopf(data, Nbus, N; backend = nothing, T = Float64, kwargs...)
    core = ExaCore(T; backend = backend)

    vr = variable(core, Nbus, N; start = ones(size(data.busarray)))
    vim = variable(core, Nbus, N;)

    pg = variable(core, size(data.gen, 1), N; lvar = repeat(data.pmin, 1, N), uvar = repeat(data.pmax, 1, N))
    qg = variable(core, size(data.gen, 1), N; lvar = repeat(data.qmin, 1, N), uvar = repeat(data.qmax, 1, N)) 

    p = variable(core, size(data.arc, 1), N; lvar = repeat(-data.rate_a, 1, N), uvar = repeat(data.rate_a, 1, N))
    q = variable(core, size(data.arc, 1), N; lvar = repeat(-data.rate_a, 1, N), uvar = repeat(data.rate_a, 1, N))

    o = objective(
        core,
        gen_cost(g, pg[g.i]) for g in data.genarray
    )

    c_ref_angle = constraint(core, c_ref_angle_rect(vr[i, t], vim[i,t]) for (i, t) in data.refarray)
 
    c_to_active_power_flow = constraint(core, c_to_active_power_flow_rect(b, p[b.f_idx, b.t],
        vr[b.f_bus, b.t], vr[b.t_bus, b.t], vim[b.f_bus, b.t], vim[b.t_bus, b.t]) for b in data.barray)
    
    c_to_reactive_power_flow = constraint(core, c_to_reactive_power_flow_rect(b, q[b.f_idx, b.t],
        vr[b.f_bus, b.t], vr[b.t_bus, b.t], vim[b.f_bus, b.t], vim[b.t_bus, b.t]) for b in data.barray)

    c_from_active_power_flow = constraint(core, c_from_active_power_flow_rect(b, p[b.t_idx, b.t],
        vr[b.f_bus, b.t], vr[b.t_bus, b.t], vim[b.f_bus, b.t], vim[b.t_bus, b.t]) for b in data.barray)

    c_from_reactive_power_flow = constraint(core, c_from_reactive_power_flow_rect(b, q[b.t_idx, b.t],
        vr[b.f_bus, b.t], vr[b.t_bus, b.t], vim[b.f_bus, b.t], vim[b.t_bus, b.t]) for b in data.barray)
    
        
    c_phase_angle_diff = constraint(
        core,
        c_phase_angle_diff_rect(b, vr[b.f_bus, b.t], vr[b.t_bus, b.t], vim[b.f_bus, b.t], vim[b.t_bus, b.t]) for b in data.barray;
        lcon = repeat(data.angmin, 1, N),
        ucon = repeat(data.angmax, 1, N),
    )

    c_active_power_balance = constraint(core, c_active_power_balance_demand_rect(b, vr[b.i, b.t], vim[b.i, b.t]) for b in data.busarray)
    c_reactive_power_balance = constraint(core, c_reactive_power_balance_demand_rect(b, vr[b.i, b.t], vim[b.i, b.t]) for b in data.busarray)
    
    c_active_power_balance_arcs = constraint!(core, c_active_power_balance, a.bus + Nbus*(a.t-1) => p[a.i, a.t] for a in data.arcarray)
    c_reactive_power_balance_arcs = constraint!(core, c_reactive_power_balance, a.bus + Nbus*(a.t-1) => q[a.i, a.t] for a in data.arcarray)

    c_active_power_balance_gen = constraint!(core, c_active_power_balance, g.bus + Nbus*(g.t-1) => -pg[g.i, g.t] for g in data.genarray)
    c_reactive_power_balance_gen = constraint!(core, c_reactive_power_balance, g.bus + Nbus*(g.t-1) => -qg[g.i, g.t] for g in data.genarray)

    c_from_thermal_limit = constraint(
        core,
        c_thermal_limit(b, p[b.f_idx, b.t], q[b.f_idx, b.t]) for b in data.barray;
        lcon = fill(-Inf, size(data.barray))
    )

    c_to_thermal_limit = constraint(
        core,
        c_thermal_limit(b, p[b.t_idx, b.t], q[b.t_idx, b.t]) for b in data.barray;
        lcon = fill(-Inf, size(data.barray))
    )

    c_voltage_magnitude = constraint(
        core, c_voltage_magnitude_rect(vr[b.i, b.t], vim[b.i, b.t])
        for b in data.busarray;
        lcon = repeat(data.vmin, 1, N).^2,
        ucon = repeat(data.vmax, 1, N).^2
    )

    c_ramp_rate = constraint(
        core,
        c_ramp(pg[g.i, g.t -1], pg[g.i, g.t]) for g in data.genarray[:, 2:N];
        lcon = repeat(-data.Δp,  1, N-1),
        ucon = repeat( data.Δp, 1, N-1),
    )

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
            c_voltage_magnitude = c_voltage_magnitude,
            c_ramp_rate = c_ramp_rate
        )
    model =ExaModel(core; kwargs...)

    return model, vars, cons
end

#Curve as input
function mpopf_model(
    filename, curve;
    N = length(curve),
    corrective_action_ratio = 0.1,
    backend = nothing,
    form = :polar,
    T = Float64,
    kwargs...,
)
    data = parse_mp_power_data(filename, N, corrective_action_ratio, backend, curve)
    
    Nbus = size(data.bus, 1)

    if form == :polar
        return build_polar_mpopf(data, Nbus, N, backend = backend, T = T, kwargs...)
    elseif form == :rect
        return build_rect_mpopf(data, Nbus, N, backend = backend, T = T, kwargs...)
    else
        error("Invalid coordinate symbol - valid options are :polar or :rect")
    end
end


#version w Pd and Qd
function mpopf_model(
    filename, active_power_data, reactive_power_data;
    pd = readdlm(active_power_data),
    qd = readdlm(reactive_power_data),
    N = size(pd, 2),
    corrective_action_ratio = 0.1,
    backend = nothing,
    form = :polar,
    T = Float64,
    kwargs...,
)
    data = parse_mp_power_data(filename, N, corrective_action_ratio, pd, qd, backend)
    
    Nbus = size(data.bus, 1)

    if form == :polar
        return build_polar_mpopf(data, Nbus, N, backend = backend, T = T, kwargs...)
    elseif form == :rect
        return build_rect_mpopf(data, Nbus, N, backend = backend, T = T, kwargs...)
    else
        error("Invalid coordinate symbol - valid options are :polar or :rect")
    end
end