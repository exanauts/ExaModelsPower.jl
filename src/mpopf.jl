using DelimitedFiles

function parse_mp_power_data(filename, N, corrective_action_ratio)

    data, dicts = parse_ac_power_data(filename)

    nbus = length(data.bus)

    empty_stor = Vector{NamedTuple{(:c, :Einit, :etac, :etad, :Srating, :Zr, :Zim, :Pexts, :Qexts, :bus, :t), Tuple{Int64, Float32, Float32, Float32, Float32, Float32, Float32, Float32, Float32, Int64, Int64}}}()

    data = (
        ;
        data...,
        refarray = [(i,t) for i in data.ref_buses, t in 1:N],
        barray = [(;b..., t = t) for b in data.branch, t in 1:N ],
        busarray = [(;b..., t = t) for b in data.bus, t in 1:N ],
        arcarray = [(;a..., t = t) for a in data.arc, t in 1:N ],
        genarray = [(;g..., t = t) for g in data.gen, t in 1:N ],
        storarray = isempty(data.storage) ? empty_data =  empty_stor : [(;s..., t = t) for s in data.storage, t in 1:N],
        Δp = corrective_action_ratio .* (data.pmax .- data.pmin)
    )
    
    return data, dicts
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

#If no storage contraints, the "build_base_xxx_mpopf" returns the final version of the mpopf
function build_base_polar_mpopf(core, data, N, Nbus)
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

    #Storage specific variables
    pstc = variable(core, size(data.storage, 1), N; lvar = zeros(size(data.storarray)), uvar = repeat(data.pcmax, 1, N))
    pstd = variable(core, size(data.storage, 1), N; lvar = zeros(size(data.storarray)), uvar = repeat(data.pdmax, 1, N))

    pst = variable(core, size(data.storage, 1), N)
    qst = variable(core, size(data.storage, 1), N)

    I2 = variable(core, size(data.storage, 1), N; lvar = zeros(size(data.storarray)))

    qint = variable(core, size(data.storage, 1), N; lvar = -repeat(data.srating, 1, N), uvar = repeat(data.srating, 1, N))

    E = variable(core, size(data.storage, 1), N; lvar = zeros(size(data.storarray)), uvar = repeat(data.emax, 1, N))

    o = objective(core, gen_cost(g, pg[g.i,g.t]) for g in data.genarray)

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
    c_reactive_power_balance_gen = constraint!(core, c_reactive_power_balance, g.bus + Nbus*(g.t-1) => -qg[g.i, g.t] for g in data.genarray)

    c_active_power_balance_stor = constraint!(core, c_active_power_balance, s.bus + Nbus*(s.t-1) => pst[s.c, s.t] for s in data.storarray)
    c_reactive_power_balance_stor = constraint!(core, c_reactive_power_balance, s.bus + Nbus*(s.t-1) => qst[s.c, s.t] for s in data.storarray)

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

    cons =  (
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

    vars = (
            va = va,
            vm = vm,
            pg = pg,
            qg = qg,
            p = p,        
            q = q, 
            pstc = pstc,
            pstd = pstd, 
            pst = pst,
            qst = qst,
            I2 = I2,
            qint = qint,
            E = E
        )

    return vars, cons
end

function build_base_rect_mpopf(core, data, N, Nbus)
    vr = variable(core, Nbus, N; start = ones(size(data.busarray)))
    vim = variable(core, Nbus, N;)

    pg = variable(core, size(data.gen, 1), N; lvar = repeat(data.pmin, 1, N), uvar = repeat(data.pmax, 1, N))
    qg = variable(core, size(data.gen, 1), N; lvar = repeat(data.qmin, 1, N), uvar = repeat(data.qmax, 1, N)) 

    p = variable(core, size(data.arc, 1), N; lvar = repeat(-data.rate_a, 1, N), uvar = repeat(data.rate_a, 1, N))
    q = variable(core, size(data.arc, 1), N; lvar = repeat(-data.rate_a, 1, N), uvar = repeat(data.rate_a, 1, N))

    #Storage specific variables
    pstc = variable(core, size(data.storage, 1), N; lvar = zeros(size(data.storarray)), uvar = repeat(data.pcmax, 1, N))
    pstd = variable(core, size(data.storage, 1), N; lvar = zeros(size(data.storarray)), uvar = repeat(data.pdmax, 1, N))

    pst = variable(core, size(data.storage, 1), N)
    qst = variable(core, size(data.storage, 1), N)

    I2 = variable(core, size(data.storage, 1), N; lvar = zeros(size(data.storarray)))

    qint = variable(core, size(data.storage, 1), N; lvar = -repeat(data.srating, 1, N), uvar = repeat(data.srating, 1, N))

    E = variable(core, size(data.storage, 1), N; lvar = zeros(size(data.storarray)), uvar = repeat(data.emax, 1, N))


    o = objective(core, gen_cost(g, pg[g.i,g.t]) for g in data.genarray)

    c_ref_angle = constraint(core, c_ref_angle_rect(vr[i, t], vim[i, t]) for (i, t) in data.refarray)
    
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
    c_active_power_balance_arcs = constraint!(core, c_reactive_power_balance, a.bus + Nbus*(a.t-1) => q[a.i, a.t] for a in data.arcarray)

    c_active_power_balance_gen = constraint!(core, c_active_power_balance, g.bus + Nbus*(g.t-1) => -pg[g.i, g.t] for g in data.genarray)
    c_reactive_power_balance_gen = constraint!(core, c_reactive_power_balance, g.bus + Nbus*(g.t-1) => -qg[g.i, g.t] for g in data.genarray)

    c_active_power_balance_stor = constraint!(core, c_active_power_balance, s.bus + Nbus*(s.t-1) => pst[s.c, s.t] for s in data.storarray)
    c_reactive_power_balance_stor = constraint!(core, c_reactive_power_balance, s.bus + Nbus*(s.t-1) => qst[s.c, s.t] for s in data.storarray)

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
            q = q,
            pstc = pstc,
            pstd = pstd, 
            pst = pst,
            qst = qst,
            I2 = I2,
            qint = qint,
            E = E
        )

    cons = cons =  (
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

    return vars, cons
end



#cc is whether to include complimentary constraint
function build_polar_mpopf(data, Nbus, N; backend = nothing, T = Float64, cc = false, kwargs...)
    core = ExaCore(T; backend = backend)

    vars, cons = build_base_polar_mpopf(core, data, N, Nbus)

    va, vm, pg, qg, p, q, pstc, pstd, pst, qst, I2, qint, E = vars

    (c_ref_angle, 
    c_to_active_power_flow, 
    c_to_reactive_power_flow, 
    c_from_active_power_flow,
    c_from_reactive_power_flow,
    c_phase_angle_diff,
    c_active_power_balance,
    c_reactive_power_balance,
    c_from_thermal_limit,
    c_to_thermal_limit,
    c_ramp_rate) = cons

    #adding storage constraints
    c_active_storage_power = constraint(core, c_active_stor_power(s, pst[s.c, s.t], pstd[s.c, s.t], pstc[s.c, s.t], I2[s.c, s.t]) for s in data.storarray)

    c_reactive_storage_power = constraint(core, c_reactive_stor_power(s, qst[s.c, s.t], qint[s.c, s.t], I2[s.c, s.t]) for s in data.storarray)

    c_ohms = constraint(core, c_ohms_polar(pst[s.c, s.t], qst[s.c, s.t], vm[s.bus, s.t], I2[s.c, s.t]) for s in data.storarray)

    c_storage_state = constraint(core, c_stor_state(s, E[s.c, s.t], E[s.c, s.t - 1], pstc[s.c, s.t], pstd[s.c, s.t]) for s in data.storarray[:, 2:N])

    c_storage_state_init = constraint(core, c_stor_state(s, E[s.c, s.t], s.Einit, pstc[s.c, s.t], pstd[s.c, s.t]) for s in data.storarray[:, 1])

    c_storage_transfer_thermal_limit  = constraint(core, c_transfer_lim(s, pst[s.c, s.t], qst[s.c, s.t]) for s in data.storarray; lcon = lcon = fill(-Inf, size(data.storarray)))

    c_discharge_thermal_limit = constraint(core, c_discharge_lim(pstd[s.c, s.t], pstc[s.c, s.t]) for s in data.storarray; lcon = -repeat(data.srating, 1, N), ucon = repeat(data.srating, 1, N))

    #Complimentarity constraint
    if cc
        c_complementarity = constraint(core, c_comp(pstc[s.c, s.t], pstd[s.c, s.t]) for s in data.storarray)
    end
    

    cons =  (
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
        c_ramp_rate = c_ramp_rate,
        c_active_storage_power = c_active_storage_power,
        c_reactive_storage_power = c_reactive_storage_power,
        c_ohms = c_ohms,
        c_storage_state = c_storage_state,
        c_storage_state_init = c_storage_state_init,
        c_storage_transfer_thermal_limit = c_storage_transfer_thermal_limit,
        c_discharge_thermal_limit = c_discharge_thermal_limit
    )

    
    model = ExaModel(core; kwargs...)
    return model, vars, cons
end

function build_rect_mpopf(data, Nbus, N; backend = nothing, T = Float64, cc = false, kwargs...)
    core = ExaCore(T; backend = backend)

    vars, cons = build_base_rect_mpopf(core, data, N, Nbus)

    vr, vim, pg, qg, p, q, pstc, pstd, pst, qst, I2, qint, E = vars

    (c_ref_angle, 
    c_to_active_power_flow, 
    c_to_reactive_power_flow, 
    c_from_active_power_flow,
    c_from_reactive_power_flow,
    c_phase_angle_diff,
    c_active_power_balance,
    c_reactive_power_balance,
    c_from_thermal_limit,
    c_to_thermal_limit,
    c_voltage_magnitude,
    c_ramp_rate) = cons
    
    #adding storage constraints
    c_active_storage_power = constraint(core, c_active_stor_power(s, pst[s.c, s.t], pstd[s.c, s.t], pstc[s.c, s.t], I2[s.c, s.t]) for s in data.storarray)

    c_reactive_storage_power = constraint(core, c_reactive_stor_power(s, qst[s.c, s.t], qint[s.c, s.t], I2[s.c, s.t]) for s in data.storarray)

    c_ohms = constraint(core, c_ohms_rect(pst[s.c, s.t], qst[s.c, s.t], vr[s.bus, s.t], vim[s.bus, s.t], I2[s.c, s.t]) for s in data.storarray)

    c_storage_state = constraint(core, c_stor_state(s, E[s.c, s.t], E[s.c, s.t - 1], pstc[s.c, s.t], pstd[s.c, s.t]) for s in data.storarray[:, 2:N])

    c_storage_state_init = constraint(core, c_stor_state(s, E[s.c, s.t], s.Einit, pstc[s.c, s.t], pstd[s.c, s.t]) for s in data.storarray[:, 1])

    c_storage_transfer_thermal_limit  = constraint(core, c_transfer_lim(s, pst[s.c, s.t], qst[s.c, s.t]) for s in data.storarray; lcon = lcon = fill(-Inf, size(data.storarray)))

    c_discharge_thermal_limit = constraint(core, c_discharge_lim(pstd[s.c, s.t], pstc[s.c, s.t]) for s in data.storarray; lcon = -repeat(data.srating, 1, N), ucon = repeat(data.srating, 1, N))

    #Complimentarity constraint
    if cc
        c_complementarity = constraint(core, c_comp(pstc[s.c, s.t], pstd[s.c, s.t]) for s in data.storarray)
    end

    cons =  (
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
        c_ramp_rate = c_ramp_rate,
        c_active_storage_power = c_active_storage_power,
        c_reactive_storage_power = c_reactive_storage_power,
        c_ohms = c_ohms,
        c_storage_state = c_storage_state,
        c_storage_state_init = c_storage_state_init,
        c_storage_transfer_thermal_limit = c_storage_transfer_thermal_limit,
        c_discharge_thermal_limit = c_discharge_thermal_limit
    )

    model = ExaModel(core; kwargs...)
    return model, vars, cons
end

#different constraints used when a function is added to remove complementarity and make charge/discharge curve smooth
function build_polar_mpopf(data, Nbus, N, discharge_func::Function; backend = nothing, T = Float64, kwargs...)
    core = ExaCore(T; backend = backend)

    vars, cons = build_base_polar_mpopf(core, data, N, Nbus)

    va, vm, pg, qg, p, q, pstc, pstd, pst, qst, I2, qint, E = vars

    (c_ref_angle, 
    c_to_active_power_flow, 
    c_to_reactive_power_flow, 
    c_from_active_power_flow,
    c_from_reactive_power_flow,
    c_phase_angle_diff,
    c_active_power_balance,
    c_reactive_power_balance,
    c_from_thermal_limit,
    c_to_thermal_limit,
    c_ramp_rate) = cons

    c_active_storage_power = constraint(core, c_active_storage_power_smooth(s, pst[s.c, s.t], pstd[s.c, s.t], I2[s.c, s.t]) for s in data.storarray)

    c_reactive_storage_power = constraint(core, c_reactive_stor_power(s, qst[s.c, s.t], qint[s.c, s.t], I2[s.c, s.t]) for s in data.storarray)

    c_ohms = constraint(core, c_ohms_polar(pst[s.c, s.t], qst[s.c, s.t], vm[s.bus, s.t], I2[s.c, s.t]) for s in data.storarray)

    c_storage_state = constraint(core, c_storage_state_smooth(s, E[s.c, s.t], E[s.c, s.t - 1], discharge_func, pstd[s.c, s.t]) for s in data.storarray[:, 2:N])

    c_storage_state_init = constraint(core, c_storage_state_smooth(s, E[s.c, s.t], s.Einit, discharge_func, pstd[s.c, s.t]) for s in data.storarray[:, 1])

    c_storage_transfer_thermal_limit  = constraint(core, c_transfer_lim(s, pst[s.c, s.t], qst[s.c, s.t]) for s in data.storarray; lcon = lcon = fill(-Inf, size(data.storarray)))

    c_discharge_thermal_limit = constraint(core, c_discharge_limit_smooth(pstd[s.c, s.t]) for s in data.storarray; lcon = -repeat(data.srating, 1, N), ucon = repeat(data.srating, 1, N))

    cons =  (
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
        c_ramp_rate = c_ramp_rate,
        c_active_storage_power = c_active_storage_power,
        c_reactive_storage_power = c_reactive_storage_power,
        c_ohms = c_ohms,
        c_storage_state = c_storage_state,
        c_storage_state_init = c_storage_state_init,
        c_storage_transfer_thermal_limit = c_storage_transfer_thermal_limit,
        c_discharge_thermal_limit = c_discharge_thermal_limit
    )

    model = ExaModel(core; kwargs...)
    return model, vars, cons
end

function build_rect_mpopf(data, Nbus, N, discharge_func::Function; backend = nothing, T = Float64, kwargs...)

    core = ExaCore(T; backend = backend)

    vars, cons = build_base_rect_mpopf(core, data, N, Nbus)

    vr, vim, pg, qg, p, q, pstc, pstd, pst, qst, I2, qint, E = vars

    (c_ref_angle, 
    c_to_active_power_flow, 
    c_to_reactive_power_flow, 
    c_from_active_power_flow,
    c_from_reactive_power_flow,
    c_phase_angle_diff,
    c_active_power_balance,
    c_reactive_power_balance,
    c_from_thermal_limit,
    c_to_thermal_limit,
    c_voltage_magnitude,
    c_ramp_rate) = cons

    c_active_storage_power = constraint(core, c_active_storage_power_smooth(s, pst[s.c, s.t], pstd[s.c, s.t], I2[s.c, s.t]) for s in data.storarray)

    c_reactive_storage_power = constraint(core, c_reactive_stor_power(s, qst[s.c, s.t], qint[s.c, s.t], I2[s.c, s.t]) for s in data.storarray)

    c_ohms = constraint(core, c_ohms_rect(pst[s.c, s.t], qst[s.c, s.t], vr[s.bus, s.t], vim[s.bus, s.t], I2[s.c, s.t]) for s in data.storarray)

    c_storage_state = constraint(core, c_storage_state_smooth(s, E[s.c, s.t], E[s.c, s.t - 1], discharge_func, pstd[s.c, s.t]) for s in data.storarray[:, 2:N])

    c_storage_state_init = constraint(core, c_storage_state_smooth(s, E[s.c, s.t], s.Einit, discharge_func, pstd[s.c, s.t]) for s in data.storarray[:, 1])

    c_storage_transfer_thermal_limit  = constraint(core, c_transfer_lim(s, pst[s.c, s.t], qst[s.c, s.t]) for s in data.storarray; lcon = lcon = fill(-Inf, size(data.storarray)))

    c_discharge_thermal_limit = constraint(core, c_discharge_limit_smooth(pstd[s.c, s.t]) for s in data.storarray; lcon = -repeat(data.srating, 1, N), ucon = repeat(data.srating, 1, N))

    cons =  (
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
        c_ramp_rate = c_ramp_rate,
        c_active_storage_power = c_active_storage_power,
        c_reactive_storage_power = c_reactive_storage_power,
        c_ohms = c_ohms,
        c_storage_state = c_storage_state,
        c_storage_state_init = c_storage_state_init,
        c_storage_transfer_thermal_limit = c_storage_transfer_thermal_limit,
        c_discharge_thermal_limit = c_discharge_thermal_limit
    )

    model = ExaModel(core; kwargs...)
    return model, vars, cons
end

function mpopf_model(
    filename, curve;
    N = length(curve),
    corrective_action_ratio = 0.1,
    backend = nothing,
    form = :polar,
    T = Float64,
    cc = false,
    kwargs...,
)

    @assert length(curve) > 0
    data, dicts = parse_mp_power_data(filename, N, corrective_action_ratio)
    update_load_data(data.busarray, curve)
    data = convert_data(data,backend)
    Nbus = size(data.bus, 1)

    if form == :polar
        return build_polar_mpopf(data, Nbus, N, backend = backend, T = T, cc = cc, kwargs...)
    elseif form == :rect
        return build_rect_mpopf(data, Nbus, N, backend = backend, T = T, cc = cc, kwargs...)
    else
        error("Invalid coordinate symbol - valid options are :polar or :rect")
    end
end

function mpopf_model(
    filename, active_power_data, reactive_power_data;
    pd = readdlm(active_power_data),
    qd = readdlm(reactive_power_data),
    N = size(pd, 2),
    corrective_action_ratio = 0.1,
    backend = nothing,
    form = :polar,
    T = Float64,
    cc = false,
    kwargs...,
)
    
    data, dicts = parse_mp_power_data(filename, N, corrective_action_ratio)
    update_load_data(data.busarray, pd, qd, data.baseMVA[], dicts.bus)
    data = convert_data(data,backend)
    Nbus = size(data.bus, 1)
    @assert Nbus == size(pd, 1)

    if form == :polar
        return build_polar_mpopf(data, Nbus, N, backend = backend, T = T, cc = cc, kwargs...)
    elseif form == :rect
        return build_rect_mpopf(data, Nbus, N, backend = backend, T = T, cc = cc, kwargs...)
    else
        error("Invalid coordinate symbol - valid options are :polar or :rect")
    end
end

#Input to discharge_func should be discharge rate (or negative charge), output should be loss in battery level
function mpopf_model(
    filename, curve, discharge_func::Function;
    N = length(curve),
    corrective_action_ratio = 0.1,
    backend = nothing,
    form = :polar,
    T = Float64,
    kwargs...,
)

    @assert length(curve) > 0
    data, dicts = parse_mp_power_data(filename, N, corrective_action_ratio)
    update_load_data(data.busarray, curve)
    data = convert_data(data,backend)
    Nbus = size(data.bus, 1)

    if form == :polar
        return build_polar_mpopf(data, Nbus, N, discharge_func, backend = backend, T = T, kwargs...)
    elseif form == :rect
        return build_rect_mpopf(data, Nbus, N, discharge_func, backend = backend, T = T, kwargs...)
    else
        error("Invalid coordinate symbol - valid options are :polar or :rect")
    end
end

function mpopf_model(
    filename, active_power_data, reactive_power_data, discharge_func::Function;
    pd = readdlm(active_power_data),
    qd = readdlm(reactive_power_data),
    N = size(pd, 2),
    corrective_action_ratio = 0.1,
    backend = nothing,
    form = :polar,
    T = Float64,
    cc = false,
    kwargs...,
)
    
    
    data, dicts = parse_mp_power_data(filename, N, corrective_action_ratio)
    update_load_data(data.busarray, pd, qd, data.baseMVA[], dicts.bus)
    data = convert_data(data,backend)
    Nbus = size(data.bus, 1)
    @assert Nbus == size(pd, 1)

    if form == :polar
        return build_polar_mpopf(data, Nbus, N, discharge_func, backend = backend, T = T, kwargs...)
    elseif form == :rect
        return build_rect_mpopf(data, Nbus, N, discharge_func, backend = backend, T = T, kwargs...)
    else
        error("Invalid coordinate symbol - valid options are :polar or :rect")
    end
end
