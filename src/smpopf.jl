using DelimitedFiles

#Curve as input
function parse_smp_power_data(filename, N, corrective_action_ratio, backend, curve)

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
        storarray = [(;s..., t = t) for s in data.storage, t in 1:N ],
        Δp = corrective_action_ratio .* (data.pmax .- data.pmin)
    )
    
    update_load_data(data.busarray, curve)
    return convert_data(data,backend)
end

#input is discharge rate (or negative charge), output is loss in battery level
function example_curve(d, srating)
    return d + .2/srating*d^2
end

function smpopf_model(
    filename, curve;
    N = length(curve),
    corrective_action_ratio = 0.1,
    backend = nothing,
    symbol = "polar",
    T = Float64,
    kwargs...,
)

    data = parse_smp_power_data(filename, N, corrective_action_ratio, backend, curve)
    Nbus = size(data.bus, 1)

    core = ExaModels.ExaCore(T; backend = backend)

    if symbol == "polar"
        va = ExaModels.variable(core, Nbus, N;)

        vm = ExaModels.variable(
            core,
            Nbus, N;
            start = ones(size(data.busarray)),
            lvar = repeat(data.vmin, 1, N),
            uvar = repeat(data.vmax, 1, N),
        )

    elseif symbol == "rect"
        vr = ExaModels.variable(core, Nbus, N; start = ones(size(data.busarray)))
        vim = ExaModels.variable(core, Nbus, N;)
    else
        error("Invalid coordinate symbol - valid options are 'polar' or 'rect'")
    end

    pg = ExaModels.variable(core, size(data.gen, 1), N; lvar = repeat(data.pmin, 1, N), uvar = repeat(data.pmax, 1, N))

    qg = ExaModels.variable(core, size(data.gen, 1), N; lvar = repeat(data.qmin, 1, N), uvar = repeat(data.qmax, 1, N)) 

    p = ExaModels.variable(core, size(data.arc, 1), N; lvar = repeat(-data.rate_a, 1, N), uvar = repeat(data.rate_a, 1, N))

    q = ExaModels.variable(core, size(data.arc, 1), N; lvar = repeat(-data.rate_a, 1, N), uvar = repeat(data.rate_a, 1, N))

    #Storage specific variables
    pstc = ExaModels.variable(core, size(data.storage, 1), N; lvar = zeros(size(data.storarray)), uvar = repeat(data.pcmax, 1, N))

    pstd = ExaModels.variable(core, size(data.storage, 1), N; lvar = zeros(size(data.storarray)), uvar = repeat(data.pdmax, 1, N))

    pst = ExaModels.variable(core, size(data.storage, 1), N)

    qst = ExaModels.variable(core, size(data.storage, 1), N)

    I2 = ExaModels.variable(core, size(data.storage, 1), N; lvar = zeros(size(data.storarray)))

    qint = ExaModels.variable(core, size(data.storage, 1), N; lvar = -repeat(data.srating, 1, N), uvar = repeat(data.srating, 1, N))

    E = ExaModels.variable(core, size(data.storage, 1), N; lvar = zeros(size(data.storarray)), uvar = repeat(data.emax, 1, N))


    o = ExaModels.objective(
        core,
        g.cost1 * pg[g.i,g.t]^2 + g.cost2 * pg[g.i,g.t] + g.cost3 for g in data.genarray
    )
    if symbol == "polar"
        c1 = ExaModels.constraint(core, va[i,t] for (i,t) in data.refarray)


        
        c2 = ExaModels.constraint(core, c2_polar(b, p[b.f_idx, b.t],
            vm[b.f_bus, b.t], vm[b.t_bus, b.t], va[b.f_bus, b.t], va[b.t_bus, b.t]) for b in data.barray)
        
        c3 = ExaModels.constraint(core, c3_polar(b, q[b.f_idx, b.t],
            vm[b.f_bus, b.t], vm[b.t_bus, b.t], va[b.f_bus, b.t], va[b.t_bus, b.t]) for b in data.barray)

        c4 = ExaModels.constraint(core, c4_polar(b, p[b.t_idx, b.t],
            vm[b.f_bus, b.t], vm[b.t_bus, b.t], va[b.f_bus, b.t], va[b.t_bus, b.t]) for b in data.barray)
        
        c5 = ExaModels.constraint(core, c5_polar(b, q[b.t_idx, b.t],
            vm[b.f_bus, b.t], vm[b.t_bus, b.t], va[b.f_bus, b.t], va[b.t_bus, b.t]) for b in data.barray)
        
        c6 = ExaModels.constraint(
            core,
            c6_polar(b, va[b.f_bus, b.t], va[b.t_bus, b.t]) for b in data.barray;
            lcon = repeat(data.angmin, 1, N),
            ucon = repeat(data.angmax, 1, N),
        )
        

        c7 = ExaModels.constraint(core, c7_polar(b, vm[b.i, b.t]) for b in data.busarray)
        c8 = ExaModels.constraint(core, c8_polar(b, vm[b.i, b.t]) for b in data.busarray)
        
        
        c9 = ExaModels.constraint(
            core,
            c9_10(b, p[b.f_idx, b.t], q[b.f_idx, b.t]) for b in data.barray;
            lcon = fill(-Inf, size(data.barray))
        )
        c10 = ExaModels.constraint(
            core,
            c9_10(b, p[b.t_idx, b.t], q[b.t_idx, b.t]) for b in data.barray;
            lcon = fill(-Inf, size(data.barray))
        )

        c15 = ExaModels.constraint(core, pst[s.c, s.t]^2 + qst[s.c, s.t]^2 - (vm[s.bus, s.t]^2)*I2[s.c, s.t] for s in data.storarray)

    elseif symbol == "rect"
        c1 = ExaModels.constraint(core, atan(vim[i, t]/vr[i, t]) for (i, t) in data.refarray)

        
        c2 = ExaModels.constraint(core, c2_rect(b, p[b.f_idx, b.t],
            vr[b.f_bus, b.t], vr[b.t_bus, b.t], vim[b.f_bus, b.t], vim[b.t_bus, b.t]) for b in data.barray)
        
        c3 = ExaModels.constraint(core, c3_rect(b, q[b.f_idx, b.t],
            vr[b.f_bus, b.t], vr[b.t_bus, b.t], vim[b.f_bus, b.t], vim[b.t_bus, b.t]) for b in data.barray)

        c4 = ExaModels.constraint(core, c4_rect(b, p[b.t_idx, b.t],
            vr[b.f_bus, b.t], vr[b.t_bus, b.t], vim[b.f_bus, b.t], vim[b.t_bus, b.t]) for b in data.barray)

        c5 = ExaModels.constraint(core, c5_rect(b, q[b.t_idx, b.t],
            vr[b.f_bus, b.t], vr[b.t_bus, b.t], vim[b.f_bus, b.t], vim[b.t_bus, b.t]) for b in data.barray)
        
            
        c6 = ExaModels.constraint(
            core,
            c6_rect(b, vr[b.f_bus, b.t], vr[b.t_bus, b.t], vim[b.f_bus, b.t], vim[b.t_bus, b.t]) for b in data.barray;
            lcon = repeat(data.angmin, 1, N),
            ucon = repeat(data.angmax, 1, N),
        )
        

        c7 = ExaModels.constraint(core, c7_rect(b, vr[b.i, b.t], vim[b.i, b.t]) for b in data.busarray)
        c8 = ExaModels.constraint(core, c8_rect(b, vr[b.i, b.t], vim[b.i, b.t]) for b in data.busarray)
        

        c9 = ExaModels.constraint(
            core,
            c9_10(b, p[b.f_idx, b.t], q[b.f_idx, b.t]) for b in data.barray;
            lcon = fill(-Inf, size(data.barray))
        )

        c10 = ExaModels.constraint(
            core,
            c9_10(b, p[b.t_idx, b.t], q[b.t_idx, b.t]) for b in data.barray;
            lcon = fill(-Inf, size(data.barray))
        )

        c11 = ExaModels.constraint(
            core, c11_rect(vr[b.i, b.t], vim[b.i, b.t])
            for b in data.busarray;
            lcon = repeat(data.vmin, 1, N).^2,
            ucon = repeat(data.vmax, 1, N).^2
        )

        c15 = ExaModels.constraint(core, pst[s.c, s.t]^2 + qst[s.c, s.t]^2 - (vr[s.bus, s.t]^2 + vim[s.bus, s.t]^2)*I2[s.c, s.t] for s in data.storarray)

    end

    c7a = ExaModels.constraint!(core, c7, a.bus + Nbus*(a.t-1) => p[a.i, a.t] for a in data.arcarray)
    c8a = ExaModels.constraint!(core, c8, a.bus + Nbus*(a.t-1) => q[a.i, a.t] for a in data.arcarray)

    c7b = ExaModels.constraint!(core, c7, g.bus + Nbus*(g.t-1) => -pg[g.i, g.t] for g in data.genarray)
    c8b = ExaModels.constraint!(core, c8, g.bus + Nbus*(g.t-1) => -qg[g.i, g.t] for g in data.genarray)

    c7c = ExaModels.constraint!(core, c7, s.bus + Nbus*(s.t-1) => pst[s.c, s.t] for s in data.storarray)
    c8c = ExaModels.constraint!(core, c8, s.bus + Nbus*(s.t-1) => qst[s.c, s.t] for s in data.storarray)


    c12 = ExaModels.constraint(
        core,
        pg[g.i, g.t -1] - pg[g.i, g.t] for g in data.genarray[:, 2:N];
        lcon = repeat(-data.Δp,  1, N-1),
        ucon = repeat( data.Δp, 1, N-1),
    )

    c13 = ExaModels.constraint(core, pst[s.c, s.t] + pstd[s.c, s.t] - pstc[s.c, s.t] - s.Pexts - s.Zr*I2[s.c, s.t] for s in data.storarray)

    c14 = ExaModels.constraint(core, qst[s.c, s.t] - qint[s.c, s.t] - s.Qexts - s.Zim*I2[s.c, s.t] for s in data.storarray)

    c16 = ExaModels.constraint(core, E[s.c, s.t] - E[s.c, s.t - 1] - (s.etac*pstc[s.c, s.t] - pstd[s.c, s.t]/s.etad) for s in data.storarray[:, 2:N])

    c17 = ExaModels.constraint(core, E[s.c, s.t] - s.Einit - (s.etac*pstc[s.c, s.t] - pstd[s.c, s.t]/s.etad) for s in data.storarray[:, 1])

    c18  = ExaModels.constraint(core, pst[s.c, s.t]^2 + qst[s.c, s.t]^2 - s.Srating^2 for s in data.storarray; lcon = lcon = fill(-Inf, size(data.storarray)))

    c19 = ExaModels.constraint(core, pstd[s.c, s.t] - pstc[s.c, s.t] for s in data.storarray; lcon = -repeat(data.srating, 1, N), ucon = repeat(data.srating, 1, N))

    #Complimentarity constraint
    c20 = ExaModels.constraint(core, pstc[s.c, s.t]*pstd[s.c, s.t] for s in data.storarray)


    if symbol == "polar"
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
    elseif symbol == "rect"
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
    end


    model =ExaModels.ExaModel(core; kwargs...)

    return model, vars
end

function smpopf_model(
    filename, curve, discharge_func::Function;
    N = length(curve),
    corrective_action_ratio = 0.1,
    backend = nothing,
    symbol = "polar",
    T = Float64,
    kwargs...,
)

    data = parse_smp_power_data(filename, N, corrective_action_ratio, backend, curve)
    Nbus = size(data.bus, 1)

    core = ExaModels.ExaCore(T; backend = backend)

    if symbol == "polar"
        va = ExaModels.variable(core, Nbus, N;)

        vm = ExaModels.variable(
            core,
            Nbus, N;
            start = ones(size(data.busarray)),
            lvar = repeat(data.vmin, 1, N),
            uvar = repeat(data.vmax, 1, N),
        )

    elseif symbol == "rect"
        vr = ExaModels.variable(core, Nbus, N; start = ones(size(data.busarray)))
        vim = ExaModels.variable(core, Nbus, N;)
    else
        error("Invalid coordinate symbol - valid options are 'polar' or 'rect'")
    end

    pg = ExaModels.variable(core, size(data.gen, 1), N; lvar = repeat(data.pmin, 1, N), uvar = repeat(data.pmax, 1, N))

    qg = ExaModels.variable(core, size(data.gen, 1), N; lvar = repeat(data.qmin, 1, N), uvar = repeat(data.qmax, 1, N)) 

    p = ExaModels.variable(core, size(data.arc, 1), N; lvar = repeat(-data.rate_a, 1, N), uvar = repeat(data.rate_a, 1, N))

    q = ExaModels.variable(core, size(data.arc, 1), N; lvar = repeat(-data.rate_a, 1, N), uvar = repeat(data.rate_a, 1, N))

    #Storage specific variables

    pstd = ExaModels.variable(core, size(data.storage, 1), N; lvar = -repeat(data.pcmax, 1, N), uvar = repeat(data.pdmax, 1, N))

    pst = ExaModels.variable(core, size(data.storage, 1), N)

    qst = ExaModels.variable(core, size(data.storage, 1), N)

    I2 = ExaModels.variable(core, size(data.storage, 1), N; lvar = zeros(size(data.storarray)))

    qint = ExaModels.variable(core, size(data.storage, 1), N; lvar = -repeat(data.srating, 1, N), uvar = repeat(data.srating, 1, N))

    E = ExaModels.variable(core, size(data.storage, 1), N; lvar = zeros(size(data.storarray)), uvar = repeat(data.emax, 1, N))

    o = ExaModels.objective(
        core,
        g.cost1 * pg[g.i,g.t]^2 + g.cost2 * pg[g.i,g.t] + g.cost3 for g in data.genarray
    )
    if symbol == "polar"
        c1 = ExaModels.constraint(core, va[i,t] for (i,t) in data.refarray)


        
        c2 = ExaModels.constraint(core, c2_polar(b, p[b.f_idx, b.t],
            vm[b.f_bus, b.t], vm[b.t_bus, b.t], va[b.f_bus, b.t], va[b.t_bus, b.t]) for b in data.barray)
        
        c3 = ExaModels.constraint(core, c3_polar(b, q[b.f_idx, b.t],
            vm[b.f_bus, b.t], vm[b.t_bus, b.t], va[b.f_bus, b.t], va[b.t_bus, b.t]) for b in data.barray)

        c4 = ExaModels.constraint(core, c4_polar(b, p[b.t_idx, b.t],
            vm[b.f_bus, b.t], vm[b.t_bus, b.t], va[b.f_bus, b.t], va[b.t_bus, b.t]) for b in data.barray)
        
        c5 = ExaModels.constraint(core, c5_polar(b, q[b.t_idx, b.t],
            vm[b.f_bus, b.t], vm[b.t_bus, b.t], va[b.f_bus, b.t], va[b.t_bus, b.t]) for b in data.barray)
        
        c6 = ExaModels.constraint(
            core,
            c6_polar(b, va[b.f_bus, b.t], va[b.t_bus, b.t]) for b in data.barray;
            lcon = repeat(data.angmin, 1, N),
            ucon = repeat(data.angmax, 1, N),
        )
        

        c7 = ExaModels.constraint(core, c7_polar(b, vm[b.i, b.t]) for b in data.busarray)
        c8 = ExaModels.constraint(core, c8_polar(b, vm[b.i, b.t]) for b in data.busarray)
        
        
        c9 = ExaModels.constraint(
            core,
            c9_10(b, p[b.f_idx, b.t], q[b.f_idx, b.t]) for b in data.barray;
            lcon = fill(-Inf, size(data.barray))
        )
        c10 = ExaModels.constraint(
            core,
            c9_10(b, p[b.t_idx, b.t], q[b.t_idx, b.t]) for b in data.barray;
            lcon = fill(-Inf, size(data.barray))
        )

        c15 = ExaModels.constraint(core, pst[s.c, s.t]^2 + qst[s.c, s.t]^2 - (vm[s.bus, s.t]^2)*I2[s.c, s.t] for s in data.storarray)

    elseif symbol == "rect"
        c1 = ExaModels.constraint(core, atan(vim[i, t]/vr[i, t]) for (i, t) in data.refarray)

        
        c2 = ExaModels.constraint(core, c2_rect(b, p[b.f_idx, b.t],
            vr[b.f_bus, b.t], vr[b.t_bus, b.t], vim[b.f_bus, b.t], vim[b.t_bus, b.t]) for b in data.barray)
        
        c3 = ExaModels.constraint(core, c3_rect(b, q[b.f_idx, b.t],
            vr[b.f_bus, b.t], vr[b.t_bus, b.t], vim[b.f_bus, b.t], vim[b.t_bus, b.t]) for b in data.barray)

        c4 = ExaModels.constraint(core, c4_rect(b, p[b.t_idx, b.t],
            vr[b.f_bus, b.t], vr[b.t_bus, b.t], vim[b.f_bus, b.t], vim[b.t_bus, b.t]) for b in data.barray)

        c5 = ExaModels.constraint(core, c5_rect(b, q[b.t_idx, b.t],
            vr[b.f_bus, b.t], vr[b.t_bus, b.t], vim[b.f_bus, b.t], vim[b.t_bus, b.t]) for b in data.barray)
        
            
        c6 = ExaModels.constraint(
            core,
            c6_rect(b, vr[b.f_bus, b.t], vr[b.t_bus, b.t], vim[b.f_bus, b.t], vim[b.t_bus, b.t]) for b in data.barray;
            lcon = repeat(data.angmin, 1, N),
            ucon = repeat(data.angmax, 1, N),
        )
        

        c7 = ExaModels.constraint(core, c7_rect(b, vr[b.i, b.t], vim[b.i, b.t]) for b in data.busarray)
        c8 = ExaModels.constraint(core, c8_rect(b, vr[b.i, b.t], vim[b.i, b.t]) for b in data.busarray)
        

        c9 = ExaModels.constraint(
            core,
            c9_10(b, p[b.f_idx, b.t], q[b.f_idx, b.t]) for b in data.barray;
            lcon = fill(-Inf, size(data.barray))
        )

        c10 = ExaModels.constraint(
            core,
            c9_10(b, p[b.t_idx, b.t], q[b.t_idx, b.t]) for b in data.barray;
            lcon = fill(-Inf, size(data.barray))
        )

        c11 = ExaModels.constraint(
            core, c11_rect(vr[b.i, b.t], vim[b.i, b.t])
            for b in data.busarray;
            lcon = repeat(data.vmin, 1, N).^2,
            ucon = repeat(data.vmax, 1, N).^2
        )

        c15 = ExaModels.constraint(core, pst[s.c, s.t]^2 + qst[s.c, s.t]^2 - (vr[s.bus, s.t]^2 + vim[s.bus, s.t]^2)*I2[s.c, s.t] for s in data.storarray)

    end

    c7a = ExaModels.constraint!(core, c7, a.bus + Nbus*(a.t-1) => p[a.i, a.t] for a in data.arcarray)
    c8a = ExaModels.constraint!(core, c8, a.bus + Nbus*(a.t-1) => q[a.i, a.t] for a in data.arcarray)

    c7b = ExaModels.constraint!(core, c7, g.bus + Nbus*(g.t-1) => -pg[g.i, g.t] for g in data.genarray)
    c8b = ExaModels.constraint!(core, c8, g.bus + Nbus*(g.t-1) => -qg[g.i, g.t] for g in data.genarray)

    c7c = ExaModels.constraint!(core, c7, s.bus + Nbus*(s.t-1) => pst[s.c, s.t] for s in data.storarray)
    c8c = ExaModels.constraint!(core, c8, s.bus + Nbus*(s.t-1) => qst[s.c, s.t] for s in data.storarray)


    c12 = ExaModels.constraint(
        core,
        pg[g.i, g.t -1] - pg[g.i, g.t] for g in data.genarray[:, 2:N];
        lcon = repeat(-data.Δp,  1, N-1),
        ucon = repeat( data.Δp, 1, N-1),
    )

    c13 = ExaModels.constraint(core, pst[s.c, s.t] + pstd[s.c, s.t] - s.Pexts - s.Zr*I2[s.c, s.t] for s in data.storarray)

    c14 = ExaModels.constraint(core, qst[s.c, s.t] - qint[s.c, s.t] - s.Qexts - s.Zim*I2[s.c, s.t] for s in data.storarray)

    c16 = ExaModels.constraint(core, E[s.c, s.t] - E[s.c, s.t - 1] + discharge_func(pstd[s.c, s.t], s.Srating) for s in data.storarray[:, 2:N])

    c17 = ExaModels.constraint(core, E[s.c, s.t] - s.Einit + discharge_func(pstd[s.c, s.t], s.Srating) for s in data.storarray[:, 1])

    c18  = ExaModels.constraint(core, pst[s.c, s.t]^2 + qst[s.c, s.t]^2 - s.Srating^2 for s in data.storarray; lcon = lcon = fill(-Inf, size(data.storarray)))

    c19 = ExaModels.constraint(core, pstd[s.c, s.t] for s in data.storarray; lcon = -repeat(data.srating, 1, N), ucon = repeat(data.srating, 1, N))

    #Complimentarity constraint
    #c20 = ExaModels.constraint(core, pstc[s.c, s.t]*pstd[s.c, s.t] for s in data.storarray)


    if symbol == "polar"
        vars = (
            va = va,
            vm = vm,
            pg = pg,
            qg = qg,
            p = p,        
            q = q, 
            pstd = pstd, 
            pst = pst,
            qst = qst,
            I2 = I2,
            qint = qint,
            E = E
        )
    elseif symbol == "rect"
        vars = (
            vr = vr,
            vim = vim,
            pg = pg,
            qg = qg,
            p = p,
            q = q,
            pstd = pstd, 
            pst = pst,
            qst = qst,
            I2 = I2,
            qint = qint,
            E = E
        )
    end


    model =ExaModels.ExaModel(core; kwargs...)

    return model, vars
end
