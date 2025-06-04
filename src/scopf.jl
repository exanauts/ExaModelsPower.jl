using GOC3Benchmark, JSON

function scopf_model(
    filename, uc_filename;
    backend = nothing,
    T = Float64,
    kwargs...
        )


    uc_data = JSON.parsefile(uc_filename)
    data = GOC3Benchmark.get_data_from_file(filename)
    sc_data, lengths = parse_sc_data(data, uc_data)

    (L_J_xf, L_J_ln, L_J_ac, L_J_dc, L_J_br, L_J_cs,
    L_J_pr, L_J_cspr, L_J_sh) = lengths

    L_T = length(sc_data.dt)
    I = length(sc_data.bus)
    Np = length(sc_data.active_reserve)
    Nq = length(sc_data.reactive_reserve)
    L_W_en_min_pr = length(sc_data.W_en_min_pr)
    L_W_en_min_cs = length(sc_data.W_en_min_cs)
    L_W_en_max_pr = length(sc_data.W_en_max_pr)
    L_W_en_max_cs = length(sc_data.W_en_max_cs)
    c_p = data.violation_cost["p_bus_vio_cost"]
    c_q = data.violation_cost["q_bus_vio_cost"]
    c_s = data.violation_cost["s_vio_cost"]
    c_e = data.violation_cost["e_vio_cost"]

    v_lvar = repeat([b.v_min for b in sc_data.bus], 1, L_T)
    v_uvar = repeat([b.v_max for b in sc_data.bus], 1, L_T)
    p_jtm_pr_uvar = [b.p_max for b in sc_data.p_jtm_flattened_pr]
    p_jtm_cs_uvar = [b.p_max for b in sc_data.p_jtm_flattened_cs]
    p_jt_fr_dc_lvar = [-dc.pdc_max for dc in sc_data.dclinearray]
    p_jt_fr_dc_uvar = [dc.pdc_max for dc in sc_data.dclinearray]
    p_jt_to_dc_lvar = [-dc.pdc_max for dc in sc_data.dclinearray]
    p_jt_to_dc_uvar = [dc.pdc_max for dc in sc_data.dclinearray]
    q_jt_fr_dc_lvar = [dc.qdc_fr_min for dc in sc_data.dclinearray]
    q_jt_fr_dc_uvar = [dc.qdc_fr_max for dc in sc_data.dclinearray]
    q_jt_to_dc_lvar = [dc.qdc_to_min for dc in sc_data.dclinearray]
    q_jt_to_dc_uvar = [dc.qdc_to_max for dc in sc_data.dclinearray]

    z_jt_on_pr = sum([pr.dt*pr.c_on*pr.u_on for pr in sc_data.prarray])
    z_jt_on_cs = sum([cs.dt*cs.c_on*cs.u_on for cs in sc_data.csarray])
    z_jt_su_pr = sum([pr.c_su*pr.u_su for pr in sc_data.prarray])
    z_jt_su_cs = sum([cs.c_su*cs.u_su for cs in sc_data.csarray])
    z_jt_su_ln = sum([ln.c_su*ln.u_su for ln in sc_data.aclbrancharray])
    z_jt_su_xf = sum([xf.c_su*xf.u_su for xf in sc_data.acxbrancharray])
    z_jt_sd_pr = sum([pr.c_sd*pr.u_sd for pr in sc_data.prarray])
    z_jt_sd_cs = sum([cs.c_sd*cs.u_sd for cs in sc_data.csarray])
    z_jt_sd_ln = sum([ln.c_sd*ln.u_sd for ln in sc_data.aclbrancharray])
    z_jt_sd_xf = sum([xf.c_sd*xf.u_sd for xf in sc_data.acxbrancharray])

    z_jt_sus_pr = 0
    for pr in sc_data.prarray
        if pr.u_su == 0
            continue
        end
        permitted_c_sus = []
        for f_ind in 1:length(pr.sus)
            if sum([tf.u_on for tf in sc_data.T_sus_jft if tf.j == pr.j && tf.t == pr.t && f_ind == tf.f], init=0) == 0
                continue
            end
            !push(permitted_c_sus, pr.sus[1])
        end
        
        min_c_sus = minimum(permitted_c_sus, init=0)
        if min_c_sus <= 0
            z_jt_sus_pr += min_c_sus
        end
    end

    z_jt_sus_cs = 0
    for cs in sc_data.csarray
        if cs.u_su == 0
            continue
        end
        permitted_c_sus = []
        for f_ind in 1:length(cs.sus)
            if sum([tf.u_on for tf in sc_data.T_sus_jft if tf.j == cs.j && tf.t == cs.t && f_ind == tf.f], init=0) == 0
                continue
            end
            !push(permitted_c_sus, cs.sus[1])
        end
        
        min_c_sus = minimum(permitted_c_sus, init=0)
        if min_c_sus <= 0
            z_jt_sus_cs += min_c_sus
        end
    end
    
    sc_data = convert_data(sc_data, backend)

    
    core = ExaCore(T; backend =backend)

    #variables are indexed j,t,k or j,t (t always second if present)

    b_jt_sh = variable(core, L_J_sh, L_T;)
    g_jt_sh = variable(core, L_J_sh, L_T;)
    #Split e_w_plus into separate sets for W_en_min and W_en_max ad for pr, cs
    #Boudns from 4.6.3 Maximum/minimum energy over multiple intervals (77)
    e_w_plus_min_pr = variable(core, L_W_en_min_pr; lvar = zeros(size(sc_data.W_en_min_pr)))
    e_w_plus_min_cs = variable(core, L_W_en_min_cs; lvar = zeros(size(sc_data.W_en_min_cs)))
    e_w_plus_max_pr = variable(core, L_W_en_max_pr; lvar = zeros(size(sc_data.W_en_max_pr)))
    e_w_plus_max_cs = variable(core, L_W_en_max_cs; lvar = zeros(size(sc_data.W_en_max_cs)))
    p_it = variable(core, I, L_T;)
    p_it_plus = variable(core, I, L_T;)
    #splitting p_jt and q_jt for shunts, producers, and consumers
    p_jt_sh = variable(core, L_J_sh, L_T;)
    p_jt_pr = variable(core, L_J_pr, L_T;)
    p_jt_cs = variable(core, L_J_cs, L_T;)
    q_jt_sh = variable(core, L_J_sh, L_T;)
    q_jt_pr = variable(core, L_J_pr, L_T;)
    q_jt_cs = variable(core, L_J_cs, L_T;)
    #Splitting p on, sd , su into pr and cs
    p_jt_on_pr = variable(core, L_J_pr, L_T;)
    p_jt_on_cs = variable(core, L_J_cs, L_T;)
    p_jt_su_pr = variable(core, L_J_pr, L_T;)
    p_jt_su_cs = variable(core, L_J_cs, L_T;)
    p_jt_sd_pr = variable(core, L_J_pr, L_T;)
    p_jt_sd_cs = variable(core, L_J_cs, L_T;)
    #p_jtm has been flattened and uses only one special index, k_flat
    #Bounds from 4.6.9 Energy cost and value (129)
    p_jtm_pr = variable(core, length(sc_data.p_jtm_flattened_pr); lvar = zeros(size(sc_data.p_jtm_flattened_pr)), uvar = p_jtm_pr_uvar)
    p_jtm_cs = variable(core, length(sc_data.p_jtm_flattened_cs); lvar = zeros(size(sc_data.p_jtm_flattened_cs)), uvar = p_jtm_cs_uvar)
    #to/from power split into ln, xf, and dc lines
    #Bounds from 4.8.4 DC lines (152-155)
    p_jt_fr_ln = variable(core, L_J_ln, L_T; start = ones(L_J_ln, L_T))
    p_jt_fr_xf = variable(core, L_J_xf, L_T; start = ones(L_J_xf, L_T))
    p_jt_fr_dc = variable(core, L_J_dc, L_T; lvar = p_jt_fr_dc_lvar, uvar = p_jt_fr_dc_uvar)
    p_jt_to_ln = variable(core, L_J_ln, L_T; start = ones(L_J_ln, L_T))
    p_jt_to_xf = variable(core, L_J_xf, L_T; start = ones(L_J_xf, L_T))
    p_jt_to_dc = variable(core, L_J_dc, L_T; lvar = p_jt_to_dc_lvar, uvar = p_jt_to_dc_uvar)
    q_jt_fr_ln = variable(core, L_J_ln, L_T; start = ones(L_J_ln, L_T))
    q_jt_fr_xf = variable(core, L_J_xf, L_T; start = ones(L_J_xf, L_T))
    q_jt_fr_dc = variable(core, L_J_dc, L_T; lvar = q_jt_fr_dc_lvar, uvar = q_jt_fr_dc_uvar)
    q_jt_to_ln = variable(core, L_J_ln, L_T; start = ones(L_J_ln, L_T))
    q_jt_to_xf = variable(core, L_J_xf, L_T; start = ones(L_J_xf, L_T))
    q_jt_to_dc = variable(core, L_J_dc, L_T; lvar = q_jt_to_dc_lvar, uvar = q_jt_to_dc_uvar)
    #p_jt rgu, rgd, scr, rru,on, rru,off, rrd,on, rrd,off and q_jt qru/qrd split into pr and cs
    #bounds from 4.6.4 Device reserve variable domains (80-89)
    p_jt_rgu_pr = variable(core, L_J_pr, L_T; lvar = zeros(size(sc_data.prarray)))
    p_jt_rgu_cs = variable(core, L_J_cs, L_T; lvar = zeros(size(sc_data.csarray)))
    p_jt_rgd_pr = variable(core, L_J_pr, L_T; lvar = zeros(size(sc_data.prarray)))
    p_jt_rgd_cs = variable(core, L_J_cs, L_T; lvar = zeros(size(sc_data.csarray)))
    p_jt_scr_pr = variable(core, L_J_pr, L_T; lvar = zeros(size(sc_data.prarray)))
    p_jt_scr_cs = variable(core, L_J_cs, L_T; lvar = zeros(size(sc_data.csarray)))
    p_jt_nsc_pr = variable(core, L_J_pr, L_T; lvar = zeros(size(sc_data.prarray)))
    p_jt_nsc_cs = variable(core, L_J_cs, L_T; lvar = zeros(size(sc_data.csarray)))
    p_jt_rru_on_pr = variable(core, L_J_pr, L_T; lvar = zeros(size(sc_data.prarray)))
    p_jt_rru_on_cs = variable(core, L_J_cs, L_T; lvar = zeros(size(sc_data.csarray)))
    p_jt_rru_off_pr = variable(core, L_J_pr, L_T; lvar = zeros(size(sc_data.prarray)))
    p_jt_rru_off_cs = variable(core, L_J_cs, L_T; lvar = zeros(size(sc_data.csarray)))
    p_jt_rrd_on_pr = variable(core, L_J_pr, L_T; lvar = zeros(size(sc_data.prarray)))
    p_jt_rrd_on_cs = variable(core, L_J_cs, L_T; lvar = zeros(size(sc_data.csarray)))
    p_jt_rrd_off_pr = variable(core, L_J_pr, L_T; lvar = zeros(size(sc_data.prarray)))
    p_jt_rrd_off_cs = variable(core, L_J_cs, L_T; lvar = zeros(size(sc_data.csarray)))
    q_jt_qru_pr = variable(core, L_J_pr, L_T; lvar = zeros(size(sc_data.prarray)))
    q_jt_qru_cs = variable(core, L_J_cs, L_T; lvar = zeros(size(sc_data.csarray)))
    q_jt_qrd_pr = variable(core, L_J_pr, L_T; lvar = zeros(size(sc_data.prarray)))
    q_jt_qrd_cs = variable(core, L_J_cs, L_T; lvar = zeros(size(sc_data.csarray)))

    
    p_nt_rgu_req = variable(core, Np, L_T;)
    p_nt_rgd_req = variable(core, Np, L_T;)
    p_nt_scr_req = variable(core, Np, L_T;)
    p_nt_nsc_req = variable(core, Np, L_T;)

    p_jt_pr_max = variable(core, L_T;)
    #Bounds from 4.3.1 Reserve shortfall domains (20-27)
    p_nt_rgu_plus = variable(core, Np, L_T; lvar = zeros(Np, L_T))
    p_nt_rgd_plus = variable(core, Np, L_T; lvar = zeros(Np, L_T))
    p_nt_scr_plus = variable(core, Np, L_T; lvar = zeros(Np, L_T))
    p_nt_nsc_plus = variable(core, Np, L_T; lvar = zeros(Np, L_T))
    p_nt_rru_plus = variable(core, Np, L_T; lvar = zeros(Np, L_T))
    p_nt_rrd_plus = variable(core, Np, L_T; lvar = zeros(Np, L_T))
    q_nt_qru_plus = variable(core, Nq, L_T; lvar = zeros(Nq, L_T))
    q_nt_qrd_plus = variable(core, Nq, L_T; lvar = zeros(Nq, L_T))

    p_t_sl = variable(core, L_T;)

    q_it = variable(core, I, L_T;)
    q_it_plus = variable(core, I, L_T;)

    #s_jt_plus split on ln and xf
    #Bounds from 4.8.1 AC branch flow limits and penalties (138)
    s_jt_plus_ln = variable(core, L_J_ln, L_T; lvar = zeros(L_J_ln, L_T))
    s_jt_plus_xf = variable(core, L_J_xf, L_T; lvar = zeros(L_J_xf, L_T))

    #Bounds form 4.2.4 Bus voltage (19)
    
    v_it = variable(core, I, L_T; lvar = v_lvar, uvar = v_uvar) 

    #skipping, for now, z_base, z_ctgavg, z_ctgmin, z_t_ctgavg, z_t_ctgmin, z_tk_ctg, z_t_t
    z_w_en_max_pr = variable(core, L_W_en_max_pr;)
    z_w_en_max_cs = variable(core, L_W_en_max_cs;)
    z_w_en_min_pr = variable(core, L_W_en_min_pr;)
    z_w_en_min_cs = variable(core, L_W_en_min_cs;)

    #split z_jt_en and on into pr and cs
    z_jt_en_pr = variable(core, L_J_pr, L_T;)
    z_jt_en_cs = variable(core, L_J_cs, L_T;)

    z_it_p = variable(core, I, L_T;)
    z_it_q = variable(core, I, L_T;)

    #z_jt_s split into ln and xf
    z_jt_s_ln = variable(core, L_J_ln, L_T;)
    z_jt_s_xf = variable(core, L_J_xf, L_T;)
    #z_jt rgu, rgd, scr, nsc, rru, rrd, qru, qrd split into pr and cs
    z_jt_rgu_pr = variable(core, L_J_pr, L_T;)
    z_jt_rgu_cs = variable(core, L_J_cs, L_T;)
    z_jt_rgd_pr = variable(core, L_J_pr, L_T;)
    z_jt_rgd_cs = variable(core, L_J_cs, L_T;)
    z_jt_scr_pr = variable(core, L_J_pr, L_T;)
    z_jt_scr_cs = variable(core, L_J_cs, L_T;)
    z_jt_nsc_pr = variable(core, L_J_pr, L_T;)
    z_jt_nsc_cs = variable(core, L_J_cs, L_T;)
    z_jt_rru_pr = variable(core, L_J_pr, L_T;)
    z_jt_rru_cs = variable(core, L_J_cs, L_T;)
    z_jt_rrd_pr = variable(core, L_J_pr, L_T;)
    z_jt_rrd_cs = variable(core, L_J_cs, L_T;)
    z_jt_qru_pr = variable(core, L_J_pr, L_T;)
    z_jt_qru_cs = variable(core, L_J_cs, L_T;)
    z_jt_qrd_pr = variable(core, L_J_pr, L_T;)
    z_jt_qrd_cs = variable(core, L_J_cs, L_T;)

    z_nt_rgu = variable(core, Np, L_T;)
    z_nt_rgd = variable(core, Np, L_T;)
    z_nt_scr = variable(core, Np, L_T;)
    z_nt_nsc = variable(core, Np, L_T;)
    z_nt_rru = variable(core, Np, L_T;)
    z_nt_rrd = variable(core, Np, L_T;)
    z_nt_qru = variable(core, Nq, L_T;)
    z_nt_qrd = variable(core, Nq, L_T;)

    θ_it = variable(core, I, L_T;)

    #split τjt and φjt into ln and xf
    τ_jt_ln = variable(core, L_J_ln, L_T; start = ones(L_J_ln, L_T))
    τ_jt_xf = variable(core, L_J_xf, L_T; start = ones(L_J_xf, L_T))
    φ_jt_ln = variable(core, L_J_ln, L_T;)
    φ_jt_xf = variable(core, L_J_xf, L_T;)


    #objective does not include contingencies rn
    #4.1 Market surplus objective
    #constraint (6-9)
    #All objectives are negative so that we can minimize

    #Removing all uc variables, which include z_on, z_su, z_sd, z_sus

    o6_t_pr = objective(core, -(-z_jt_en_pr[pr.j_pr, pr.t] 
                        - (z_jt_rgu_pr[pr.j_pr, pr.t] + z_jt_rgd_pr[pr.j_pr, pr.t] + z_jt_scr_pr[pr.j_pr, pr.t] + z_jt_nsc_pr[pr.j_pr, pr.t] + z_jt_rru_pr[pr.j_pr, pr.t] + 
                        z_jt_rrd_pr[pr.j_pr, pr.t] + z_jt_qru_pr[pr.j_pr, pr.t] +z_jt_qrd_pr[pr.j_pr, pr.t])) for pr in sc_data.prarray)
    o6_t_cs = objective(core, -(z_jt_en_cs[cs.j_cs, cs.t] 
                        - (z_jt_rgu_cs[cs.j_cs, cs.t] + z_jt_rgd_cs[cs.j_cs, cs.t] + z_jt_scr_cs[cs.j_cs, cs.t] + z_jt_nsc_cs[cs.j_cs, cs.t] + z_jt_rru_cs[cs.j_cs, cs.t] + 
                        z_jt_rrd_cs[cs.j_cs, cs.t] + z_jt_qru_cs[cs.j_cs, cs.t] +z_jt_qrd_cs[cs.j_cs, cs.t])) for cs in sc_data.csarray)
    o6_t_ln = objective(core, -(- z_jt_s_ln[ln.j_ln, ln.t]) for ln in sc_data.aclbrancharray)
    o6_t_xf = objective(core, -(- z_jt_s_xf[xf.j_xf, xf.t]) for xf in sc_data.acxbrancharray)
    o6_t_i = objective(core, -(-(z_it_p[b.i, b.t] + z_it_q[b.i, b.t])) for b in sc_data.busarray)
    o6_t_Np = objective(core, -(-(z_nt_rgu[n.n_p, n.t] + z_nt_rgd[n.n_p, n.t] + z_nt_scr[n.n_p, n.t] + z_nt_nsc[n.n_p, n.t] + z_nt_rru[n.n_p, n.t] + z_nt_rrd[n.n_p, n.t])) for n in sc_data.preservearray)
    o6_t_Nq = objective(core, -(-(z_nt_qru[n.n_q, n.t] + z_nt_qrd[n.n_q, n.t])) for n in sc_data.qreservearray)
    o6_en_max_pr = objective(core, z_w_en_max_pr[w] for w in 1:L_W_en_max_pr)
    o6_en_max_cs = objective(core, z_w_en_max_cs[w] for w in 1:L_W_en_max_cs)
    o6_en_min_pr = objective(core, z_w_en_min_pr[w] for w in 1:L_W_en_min_pr)
    o6_en_min_cs = objective(core, z_w_en_min_cs[w] for w in 1:L_W_en_min_cs)


    #4.2.1 Bus power mismatch and penalized mismatch definitions
    c_11 = constraint(core, p_it_plus[b.i, b.t] - p_it[b.i, b.t] for b in sc_data.busarray; ucon = fill(Inf, size(sc_data.busarray)))
    c_12 = constraint(core, p_it_plus[b.i, b.t] + p_it[b.i, b.t] for b in sc_data.busarray; ucon = fill(Inf, size(sc_data.busarray)))
    c_13 = constraint(core, q_it_plus[b.i, b.t] - q_it[b.i, b.t] for b in sc_data.busarray; ucon = fill(Inf, size(sc_data.busarray)))
    c_14 = constraint(core, q_it_plus[b.i, b.t] + q_it[b.i, b.t] for b in sc_data.busarray; ucon = fill(Inf, size(sc_data.busarray)))
    #4.2.2 Bus pwoer mismatch penalty
    c15 = constraint(core, z_it_p[b.i, b.t] - b.dt*c_p*p_it_plus[b.i, b.t] for b in sc_data.busarray)
    c16 = constraint(core, z_it_q[b.i, b.t] - b.dt*c_q*q_it_plus[b.i, b.t] for b in sc_data.busarray)
    #4.2.3 Bus real and reactive power balance
    c17 = constraint(core, p_it[b.i, b.t] for b in sc_data.busarray)
    c17_pr = constraint!(core, c17, pr.bus + I*(pr.t-1) => p_jt_pr[pr.j_pr, pr.t] for pr in sc_data.prarray)
    c17_cs = constraint!(core, c17, cs.bus + I*(cs.t-1) => -p_jt_cs[cs.j_cs, cs.t] for cs in sc_data.csarray)
    c17_sh = constraint!(core, c17, sh.bus + I*(sh.t-1) => -p_jt_sh[sh.j_sh, sh.t] for sh in sc_data.shuntarray)
    #Reminder: fr and to split for ln, xf, and dc
    c17_fr_ln = constraint!(core, c17, ln.fr_bus + I*(ln.t-1) => -p_jt_fr_ln[ln.j_ln, ln.t] for ln in sc_data.aclbrancharray)
    c17_fr_xf = constraint!(core, c17, xf.fr_bus + I*(xf.t-1) => -p_jt_fr_xf[xf.j_xf, xf.t] for xf in sc_data.acxbrancharray)
    c17_fr_dc = constraint!(core, c17, dc.fr_bus + I*(dc.t-1) => -p_jt_fr_dc[dc.j_dc, dc.t] for dc in sc_data.dclinearray)
    c17_to_ln = constraint!(core, c17, ln.to_bus + I*(ln.t-1) => -p_jt_to_ln[ln.j_ln, ln.t] for ln in sc_data.aclbrancharray)
    c17_to_xf = constraint!(core, c17, xf.to_bus + I*(xf.t-1) => -p_jt_to_xf[xf.j_xf, xf.t] for xf in sc_data.acxbrancharray)
    c17_to_dc = constraint!(core, c17, dc.to_bus + I*(dc.t-1) => -p_jt_to_dc[dc.j_dc, dc.t] for dc in sc_data.dclinearray)
    c18 = constraint(core, q_it[b.i, b.t] for b in sc_data.busarray)
    c18_pr = constraint!(core, c18, pr.bus + I*(pr.t-1) => q_jt_pr[pr.j_pr, pr.t] for pr in sc_data.prarray)
    c18_cs = constraint!(core, c18, cs.bus + I*(cs.t-1) => -q_jt_cs[cs.j_cs, cs.t] for cs in sc_data.csarray)
    c18_sh = constraint!(core, c18, sh.bus + I*(sh.t-1) => -q_jt_sh[sh.j_sh, sh.t] for sh in sc_data.shuntarray)
    #Reminder: fr and to split for ln, xf, and dc
    c18_fr_ln = constraint!(core, c18, ln.fr_bus + I*(ln.t-1) => -q_jt_fr_ln[ln.j_ln, ln.t] for ln in sc_data.aclbrancharray)
    c18_fr_xf = constraint!(core, c18, xf.fr_bus + I*(xf.t-1) => -q_jt_fr_xf[xf.j_xf, xf.t] for xf in sc_data.acxbrancharray)
    c18_fr_dc = constraint!(core, c18, dc.fr_bus + I*(dc.t-1) => -q_jt_fr_dc[dc.j_dc, dc.t] for dc in sc_data.dclinearray)
    c18_to_ln = constraint!(core, c18, ln.to_bus + I*(ln.t-1) => -q_jt_to_ln[ln.j_ln, ln.t] for ln in sc_data.aclbrancharray)
    c18_to_xf = constraint!(core, c18, xf.to_bus + I*(xf.t-1) => -q_jt_to_xf[xf.j_xf, xf.t] for xf in sc_data.acxbrancharray)
    c18_to_dc = constraint!(core, c18, dc.to_bus + I*(dc.t-1) => -q_jt_to_dc[dc.j_dc, dc.t] for dc in sc_data.dclinearray)

    #4.3.2 Reserve shortfall penalties
    c28 = constraint(core, z_nt_rgu[n.n_p, n.t] - n.dt*n.c_rgu*p_nt_rgu_plus[n.n_p, n.t] for n in sc_data.preservearray)
    c29 = constraint(core, z_nt_rgd[n.n_p, n.t] - n.dt*n.c_rgd*p_nt_rgd_plus[n.n_p, n.t] for n in sc_data.preservearray)
    c30 = constraint(core, z_nt_scr[n.n_p, n.t] - n.dt*n.c_scr*p_nt_scr_plus[n.n_p, n.t] for n in sc_data.preservearray)
    c31 = constraint(core, z_nt_nsc[n.n_p, n.t] - n.dt*n.c_nsc*p_nt_nsc_plus[n.n_p, n.t] for n in sc_data.preservearray)
    c32 = constraint(core, z_nt_rru[n.n_p, n.t] - n.dt*n.c_rru*p_nt_rru_plus[n.n_p, n.t] for n in sc_data.preservearray)
    c33 = constraint(core, z_nt_rrd[n.n_p, n.t] - n.dt*n.c_rrd*p_nt_rrd_plus[n.n_p, n.t] for n in sc_data.preservearray)
    c34 = constraint(core, z_nt_qru[n.n_q, n.t] - n.dt*n.c_qru*q_nt_qru_plus[n.n_q, n.t] for n in sc_data.qreservearray)
    c35 = constraint(core, z_nt_qrd[n.n_q, n.t] - n.dt*n.c_qrd*q_nt_qrd_plus[n.n_q, n.t] for n in sc_data.qreservearray)

    #4.3.3 Reserve requirements
    c36 = constraint(core, p_nt_rgu_req[n.n_p, n.t]/n.σ_rgu for n in sc_data.preservearray)
    c36_cs = constraint!(core, c36, cs.n + Np*(cs.t-1) => -p_jt_cs[cs.j_cs, cs.t] for cs in sc_data.preservesetarray_cs)
    c37 = constraint(core, p_nt_rgd_req[n.n_p, n.t]/n.σ_rgd for n in sc_data.preservearray)
    c37_cs = constraint!(core, c37, cs.n + Np*(cs.t-1) => -p_jt_cs[cs.j_cs, cs.t] for cs in sc_data.preservesetarray_cs)
    #assuming c_scr and c_nsc are always positive
    cmax38 = constraint(core, p_jt_pr_max[pr.t] - p_jt_pr[pr.j_pr, pr.t] for pr in sc_data.prarray; ucon = fill(Inf, size(sc_data.prarray)))
    c38 = constraint(core, p_nt_scr_req[n.n_p, n.t] - n.σ_scr*p_jt_pr_max[n.t] for n in sc_data.preservearray)
    c39 = constraint(core, p_nt_nsc_req[n.n_p, n.t] - n.σ_nsc*p_jt_pr_max[n.t] for n in sc_data.preservearray)

    #4.3.4 Reserve balance
    #Reminder, p and q sets have been split up for pr and cs
    c40 = constraint(core, p_nt_rgu_plus[n.n_p, n.t] - p_nt_rgu_req[n.n_p, n.t] for n in sc_data.preservearray; ucon = fill(Inf, size(sc_data.preservearray)))
    c40_pr = constraint!(core, c40, pr.n_p + Np*(pr.t-1) => p_jt_rgu_pr[pr.j_pr, pr.t] for pr in sc_data.preservesetarray_pr)
    c40_cs = constraint!(core, c40, cs.n_p + Np*(cs.t-1) => p_jt_rgu_cs[cs.j_cs, cs.t] for cs in sc_data.preservesetarray_cs)
    c41 = constraint(core, p_nt_rgd_plus[n.n_p, n.t] - p_nt_rgd_req[n.n_p, n.t] for n in sc_data.preservearray; ucon = fill(Inf, size(sc_data.preservearray)))
    c41_pr = constraint!(core, c41, pr.n_p + Np*(pr.t-1) => p_jt_rgd_pr[pr.j_pr, pr.t] for pr in sc_data.preservesetarray_pr)
    c41_cs = constraint!(core, c41, cs.n_p + Np*(cs.t-1) => p_jt_rgd_cs[cs.j_cs, cs.t] for cs in sc_data.preservesetarray_cs)
    c42 = constraint(core, p_nt_scr_plus[n.n_p, n.t] - p_nt_rgu_req[n.n_p, n.t] - p_nt_scr_req[n.n_p, n.t] for n in sc_data.preservearray; ucon = fill(Inf, size(sc_data.preservearray)))
    c42_pr = constraint!(core, c42, pr.n_p + Np*(pr.t-1) => p_jt_rgu_pr[pr.j_pr, pr.t] + p_jt_scr_pr[pr.j_pr, pr.t] for pr in sc_data.preservesetarray_pr)
    c42_cs = constraint!(core, c42, cs.n_p + Np*(cs.t-1) => p_jt_rgu_cs[cs.j_cs, cs.t] + p_jt_scr_cs[cs.j_cs, cs.t] for cs in sc_data.preservesetarray_cs)
    c43 = constraint(core, p_nt_nsc_plus[n.n_p, n.t] - p_nt_rgu_req[n.n_p, n.t] - p_nt_scr_req[n.n_p, n.t] - p_nt_nsc_req[n.n_p, n.t] for n in sc_data.preservearray; ucon = fill(Inf, size(sc_data.preservearray)))
    c43_pr = constraint!(core, c43, pr.n_p + Np*(pr.t-1) => p_jt_rgu_pr[pr.j_pr, pr.t] + p_jt_scr_pr[pr.j_pr, pr.t] + p_jt_nsc_pr[pr.j_pr, pr.t] for pr in sc_data.preservesetarray_pr)
    c43_cs = constraint!(core, c43, cs.n_p + Np*(cs.t-1) => p_jt_rgu_cs[cs.j_cs, cs.t] + p_jt_scr_cs[cs.j_cs, cs.t] + p_jt_nsc_cs[cs.j_cs, cs.t] for cs in sc_data.preservesetarray_cs)
    c44 = constraint(core, p_nt_rru_plus[n.n_p, n.t] - n.p_rru_min for n in sc_data.preservearray; ucon = fill(Inf, size(sc_data.preservearray)))
    c44_pr = constraint!(core, c44, pr.n_p + Np*(pr.t-1) => p_jt_rru_on_pr[pr.j_pr, pr.t] + p_jt_rru_off_pr[pr.j_pr, pr.t] for pr in sc_data.preservesetarray_pr)
    c44_cs = constraint!(core, c44, cs.n_p + Np*(cs.t-1) => p_jt_rru_on_cs[cs.j_cs, cs.t] + p_jt_rru_off_cs[cs.j_cs, cs.t] for cs in sc_data.preservesetarray_cs)
    c45 = constraint(core, p_nt_rrd_plus[n.n_p, n.t] - n.p_rrd_min for n in sc_data.preservearray; ucon = fill(Inf, size(sc_data.preservearray)))
    c45_pr = constraint!(core, c45, pr.n_p + Np*(pr.t-1) => p_jt_rrd_on_pr[pr.j_pr, pr.t] + p_jt_rrd_off_pr[pr.j_pr, pr.t] for pr in sc_data.preservesetarray_pr)
    c45_cs = constraint!(core, c45, cs.n_p + Np*(cs.t-1) => p_jt_rrd_on_cs[cs.j_cs, cs.t] + p_jt_rrd_off_cs[cs.j_cs, cs.t] for cs in sc_data.preservesetarray_cs)
    c46 = constraint(core, q_nt_qru_plus[n.n_q, n.t] - n.q_qru_min for n in sc_data.qreservearray; ucon = fill(Inf, size(sc_data.qreservearray)))
    c46_pr = constraint!(core, c46, pr.n_q + Nq*(pr.t-1) => q_jt_qru_pr[pr.j_pr, pr.t] for pr in sc_data.qreservesetarray_pr)
    c46_cs = constraint!(core, c46, cs.n_q + Nq*(cs.t-1) => q_jt_qru_cs[cs.j_cs, cs.t] for cs in sc_data.qreservesetarray_cs)
    c47 = constraint(core, q_nt_qrd_plus[n.n_q, n.t] - n.q_qrd_min for n in sc_data.qreservearray; ucon = fill(Inf, size(sc_data.qreservearray)))
    c47_pr = constraint!(core, c47, pr.n_q + Nq*(pr.t-1) => q_jt_qrd_pr[pr.j_pr, pr.t] for pr in sc_data.qreservesetarray_pr)
    c47_cs = constraint!(core, c47, cs.n_q + Nq*(cs.t-1) => q_jt_qrd_cs[cs.j_cs, cs.t] for cs in sc_data.qreservesetarray_cs)

    #skipping constraints 48-67. Assume unit commitment is always satisfied
    
    #4.6.1 Producing and consuming device startup, shutdown, and dispatchable power
    #p_jt variants are split on pr and cs
    c68_pr = constraint(core, p_jt_on_pr[pr.j_pr, pr.t] + p_jt_su_pr[pr.j_pr, pr.t] + p_jt_sd_pr[pr.j_pr, pr.t] - p_jt_pr[pr.j_pr, pr.t] for pr in sc_data.prarray)
    c68_cs = constraint(core, p_jt_on_cs[cs.j_cs, cs.t] + p_jt_su_cs[cs.j_cs, cs.t] + p_jt_sd_cs[cs.j_cs, cs.t] - p_jt_cs[cs.j_cs, cs.t] for cs in sc_data.csarray)
    c69_pr = constraint(core, p_jt_su_pr[pr.j_pr, pr.t] - pr.sum_T_supc_pr_jt for pr in sc_data.prarray)
    c69_cs = constraint(core, p_jt_su_cs[cs.j_cs, cs.t] - cs.sum_T_supc_cs_jt for cs in sc_data.csarray)
    c70_pr = constraint(core, p_jt_sd_pr[pr.j_pr, pr.t] - pr.sum_T_sdpc_pr_jt for pr in sc_data.prarray)
    c70_cs = constraint(core, p_jt_sd_cs[cs.j_cs, cs.t] - cs.sum_T_sdpc_cs_jt for cs in sc_data.csarray)

    #4.6.2 Ramping limits
    #p split for pr and cs
    c71_pr = constraint(core, p_jt_pr[pr.j_pr, pr.t] - pr.p_0 - pr.dt*(pr.p_ru*(pr.u_on-pr.u_su) + pr.p_ru_su*(pr.u_su+1-pr.u_on)) for pr in sc_data.prarray[1:L_J_pr];
    lcon = fill(-Inf, size(sc_data.prarray[1:L_J_pr])))
    c71_cs = constraint(core, p_jt_cs[cs.j_cs, cs.t] - cs.p_0 - cs.dt*(cs.p_ru*(cs.u_on-cs.u_su) + cs.p_ru_su*(cs.u_su+1-cs.u_on)) for cs in sc_data.csarray[1:L_J_cs];
    lcon  = fill(-Inf, size(sc_data.csarray[1:L_J_cs])))
    c72_pr = constraint(core, p_jt_pr[pr.j_pr, pr.t] - p_jt_pr[pr.j_pr, pr.t-1] - pr.dt*(pr.p_ru*(pr.u_on-pr.u_su) + pr.p_ru_su*(pr.u_su+1-pr.u_on)) for pr in sc_data.prarray[L_J_pr+1:end];
    lcon = fill(-Inf, size(sc_data.prarray[L_J_pr+1:end])))
    c72_cs = constraint(core, p_jt_cs[cs.j_cs, cs.t] - p_jt_cs[cs.j_cs, cs.t-1] - cs.dt*(cs.p_ru*(cs.u_on-cs.u_su) + cs.p_ru_su*(cs.u_su+1-cs.u_on)) for cs in sc_data.csarray[L_J_cs+1:end];
    lcon = fill(-Inf, size(sc_data.csarray[L_J_cs+1:end])))
    c73_pr = constraint(core, p_jt_pr[pr.j_pr, pr.t] - pr.p_0 + pr.dt*(pr.p_rd*pr.u_on+pr.p_rd_sd*(1-pr.u_on)) for pr in sc_data.prarray[1:L_J_pr];
    ucon = fill(Inf, size(sc_data.prarray[1:L_J_pr])))
    c73_cs = constraint(core, p_jt_cs[cs.j_cs, cs.t] - cs.p_0 + cs.dt*(cs.p_rd*cs.u_on+cs.p_rd_sd*(1-cs.u_on)) for cs in sc_data.csarray[1:L_J_cs];
    ucon = fill(Inf, size(sc_data.csarray[1:L_J_cs])))
    c74_pr = constraint(core, p_jt_pr[pr.j_pr, pr.t] - p_jt_pr[pr.j_pr, pr.t-1] + pr.dt*(pr.p_rd*pr.u_on+pr.p_rd_sd*(1-pr.u_on)) for pr in sc_data.prarray[L_J_pr+1:end];
    ucon = fill(Inf, size(sc_data.prarray[L_J_pr+1:end])))
    c74_cs = constraint(core, p_jt_cs[cs.j_cs, cs.t] - p_jt_cs[cs.j_cs, cs.t-1] + cs.dt*(cs.p_rd*cs.u_on+cs.p_rd_sd*(1-cs.u_on)) for cs in sc_data.csarray[L_J_cs+1:end];
    ucon = fill(Inf, size(sc_data.csarray[L_J_cs+1:end])))

    #4.6.3 Maximum/minimum energy over multiple intervals
    #J_pr,cs has been split for pr and cs
    c75_pr = constraint(core, e_w_plus_max_pr[w.w_en_max_pr_ind] + w.e_max for w in sc_data.W_en_max_pr; ucon = fill(Inf, size(sc_data.W_en_max_pr)))
    c75_pr_a = constraint!(core, c75_pr, t.w_en_max_pr_ind => -t.dt*p_jt_pr[t.j_pr, t.t] for t in sc_data.T_w_en_max_pr)
    c75_cs = constraint(core, e_w_plus_max_cs[w.w_en_max_cs_ind] + w.e_max for w in sc_data.W_en_max_cs; ucon = fill(Inf, size(sc_data.W_en_max_cs)))
    c75_cs_a = constraint!(core, c75_cs, t.w_en_max_cs_ind => -t.dt*p_jt_cs[t.j_cs, t.t] for t in sc_data.T_w_en_max_cs)
    c76_pr = constraint(core, e_w_plus_min_pr[w.w_en_min_pr_ind] - w.e_min for w in sc_data.W_en_min_pr; lcon = fill(-Inf, size(sc_data.W_en_min_pr)))
    c76_pr_a = constraint!(core, c76_pr, t.w_en_min_pr_ind => -t.dt*p_jt_pr[t.j_pr, t.t] for t in sc_data.T_w_en_min_pr)
    c76_cs = constraint(core, e_w_plus_min_cs[w.w_en_min_cs_ind] - w.e_min for w in sc_data.W_en_min_cs; lcon = fill(-Inf, size(sc_data.W_en_min_cs)))
    c76_cs_a = constraint!(core, c76_cs, t.w_en_min_cs_ind => -t.dt*p_jt_cs[t.j_cs, t.t] for t in sc_data.T_w_en_min_cs)
    c78_pr = constraint(core, z_w_en_max_pr[w.w_en_max_pr_ind] - c_e*e_w_plus_max_pr[w.w_en_max_pr_ind] for w in sc_data.W_en_max_pr)
    c78_cs = constraint(core, z_w_en_max_cs[w.w_en_max_cs_ind] - c_e*e_w_plus_max_cs[w.w_en_max_cs_ind] for w in sc_data.W_en_max_cs)
    c79_pr = constraint(core, z_w_en_min_pr[w.w_en_min_pr_ind] - c_e*e_w_plus_min_pr[w.w_en_min_pr_ind] for w in sc_data.W_en_min_pr)
    c79_cs = constraint(core, z_w_en_min_cs[w.w_en_min_cs_ind] - c_e*e_w_plus_min_cs[w.w_en_min_cs_ind] for w in sc_data.W_en_min_cs)

    #4.6.5 Device reserve costs
    #p_jt split into pr and cs
    c90_pr = constraint(core, z_jt_rgu_pr[pr.j_pr, pr.t] - pr.dt*pr.c_rgu*p_jt_rgu_pr[pr.j_pr, pr.t] for pr in sc_data.prarray)
    c90_cs = constraint(core, z_jt_rgu_cs[cs.j_cs, cs.t] - cs.dt*cs.c_rgu*p_jt_rgu_cs[cs.j_cs, cs.t] for cs in sc_data.csarray)
    c91_pr = constraint(core, z_jt_rgd_pr[pr.j_pr, pr.t] - pr.dt*pr.c_rgd*p_jt_rgd_pr[pr.j_pr, pr.t] for pr in sc_data.prarray)
    c91_cs = constraint(core, z_jt_rgd_cs[cs.j_cs, cs.t] - cs.dt*cs.c_rgd*p_jt_rgd_cs[cs.j_cs, cs.t] for cs in sc_data.csarray)
    c92_pr = constraint(core, z_jt_scr_pr[pr.j_pr, pr.t] - pr.dt*pr.c_scr*p_jt_scr_pr[pr.j_pr, pr.t] for pr in sc_data.prarray)
    c92_cs = constraint(core, z_jt_scr_cs[cs.j_cs, cs.t] - cs.dt*cs.c_scr*p_jt_scr_cs[cs.j_cs, cs.t] for cs in sc_data.csarray)
    c93_pr = constraint(core, z_jt_nsc_pr[pr.j_pr, pr.t] - pr.dt*pr.c_nsc*p_jt_nsc_pr[pr.j_pr, pr.t] for pr in sc_data.prarray)
    c93_cs = constraint(core, z_jt_nsc_cs[cs.j_cs, cs.t] - cs.dt*cs.c_nsc*p_jt_nsc_cs[cs.j_cs, cs.t] for cs in sc_data.csarray)
    c94_pr = constraint(core, z_jt_rru_pr[pr.j_pr, pr.t] - pr.dt*(pr.c_rru_on*p_jt_rru_on_pr[pr.j_pr, pr.t] + pr.c_rru_off*p_jt_rru_off_pr[pr.j_pr, pr.t]) for pr in sc_data.prarray)
    c94_cs = constraint(core, z_jt_rru_cs[cs.j_cs, cs.t] - cs.dt*(cs.c_rru_on*p_jt_rru_on_cs[cs.j_cs, cs.t] + cs.c_rru_off*p_jt_rru_off_cs[cs.j_cs, cs.t]) for cs in sc_data.csarray)
    c95_pr = constraint(core, z_jt_rrd_pr[pr.j_pr, pr.t] - pr.dt*(pr.c_rrd_on*p_jt_rrd_on_pr[pr.j_pr, pr.t] + pr.c_rrd_off*p_jt_rrd_off_pr[pr.j_pr, pr.t]) for pr in sc_data.prarray)
    c95_cs = constraint(core, z_jt_rrd_cs[cs.j_cs, cs.t] - cs.dt*(cs.c_rrd_on*p_jt_rrd_on_cs[cs.j_cs, cs.t] + cs.c_rrd_off*p_jt_rrd_off_cs[cs.j_cs, cs.t]) for cs in sc_data.csarray)
    c96_pr = constraint(core, z_jt_qru_pr[pr.j_pr, pr.t] - pr.dt*pr.c_qru*q_jt_qru_pr[pr.j_pr, pr.t] for pr in sc_data.prarray)
    c96_cs = constraint(core, z_jt_qru_cs[cs.j_cs, cs.t] - cs.dt*cs.c_qru*q_jt_qru_cs[cs.j_cs, cs.t] for cs in sc_data.csarray)
    c97_pr = constraint(core, z_jt_qrd_pr[pr.j_pr, pr.t] - pr.dt*pr.c_qrd*q_jt_qrd_pr[pr.j_pr, pr.t] for pr in sc_data.prarray)
    c97_cs = constraint(core, z_jt_qrd_cs[cs.j_cs, cs.t] - cs.dt*cs.c_qrd*q_jt_qrd_cs[cs.j_cs, cs.t] for cs in sc_data.csarray)

    #4.6.6 Absolute reserve limits, based on ramp rates
    c98_pr = constraint(core, p_jt_rgu_pr[pr.j_pr, pr.t] - pr.p_rgu_max*pr.u_on for pr in sc_data.prarray; lcon = fill(-Inf, size(sc_data.prarray)))
    c98_cs = constraint(core, p_jt_rgu_cs[cs.j_cs, cs.t] - cs.p_rgu_max*cs.u_on for cs in sc_data.csarray; lcon = fill(-Inf, size(sc_data.csarray)))
    c99_pr = constraint(core, p_jt_rgd_pr[pr.j_pr, pr.t] - pr.p_rgd_max*pr.u_on for pr in sc_data.prarray; lcon = fill(-Inf, size(sc_data.prarray)))
    c99_cs = constraint(core, p_jt_rgd_cs[cs.j_cs, cs.t] - cs.p_rgd_max*cs.u_on for cs in sc_data.csarray; lcon = fill(-Inf, size(sc_data.csarray)))
    c100_pr = constraint(core, p_jt_rgu_pr[pr.j_pr, pr.t] + p_jt_scr_pr[pr.j_pr, pr.t] - pr.p_scr_max*pr.u_on for pr in sc_data.prarray; lcon = fill(-Inf, size(sc_data.prarray)))
    c100_cs = constraint(core, p_jt_rgu_cs[cs.j_cs, cs.t] + p_jt_scr_cs[cs.j_cs, cs.t] - cs.p_scr_max*cs.u_on for cs in sc_data.csarray; lcon = fill(-Inf, size(sc_data.csarray)))
    c101_pr = constraint(core, p_jt_nsc_pr[pr.j_pr, pr.t] - pr.p_nsc_max*(1-pr.u_on) for pr in sc_data.prarray; lcon = fill(-Inf, size(sc_data.prarray)))
    c101_cs = constraint(core, p_jt_nsc_cs[cs.j_cs, cs.t] - cs.p_nsc_max*(1-cs.u_on) for cs in sc_data.csarray; lcon = fill(-Inf, size(sc_data.csarray)))
    c102_pr = constraint(core, p_jt_rgu_pr[pr.j_pr, pr.t] + p_jt_scr_pr[pr.j_pr, pr.t] + p_jt_rru_on_pr[pr.j_pr, pr.t] - pr.p_rru_on_max*pr.u_on for pr in sc_data.prarray; lcon = fill(-Inf, size(sc_data.prarray)))
    c102_cs = constraint(core, p_jt_rgu_cs[cs.j_cs, cs.t] + p_jt_scr_cs[cs.j_cs, cs.t] + p_jt_rru_on_cs[cs.j_cs, cs.t] - cs.p_rru_on_max*cs.u_on for cs in sc_data.csarray; lcon = fill(-Inf, size(sc_data.csarray)))
    c103_pr = constraint(core, p_jt_nsc_pr[pr.j_pr, pr.t] + p_jt_rru_off_pr[pr.j_pr, pr.t] - pr.p_rru_off_max*(1-pr.u_on) for pr in sc_data.prarray; lcon = fill(-Inf, size(sc_data.prarray)))
    c103_cs = constraint(core, p_jt_nsc_cs[cs.j_cs, cs.t] + p_jt_rru_off_cs[cs.j_cs, cs.t] - cs.p_rru_off_max*(1-cs.u_on) for cs in sc_data.csarray; lcon = fill(-Inf, size(sc_data.csarray)))
    c104_pr = constraint(core, p_jt_rgd_pr[pr.j_pr, pr.t] + p_jt_rrd_on_pr[pr.j_pr, pr.t] - pr.p_rrd_on_max*pr.u_on for pr in sc_data.prarray; lcon = fill(-Inf, size(sc_data.prarray)))
    c104_cs = constraint(core, p_jt_rgd_cs[cs.j_cs, cs.t] + p_jt_rrd_on_cs[cs.j_cs, cs.t] - cs.p_rrd_on_max*cs.u_on for cs in sc_data.csarray; lcon = fill(-Inf, size(sc_data.csarray)))
    c105_pr = constraint(core, p_jt_rrd_off_pr[pr.j_pr, pr.t] - pr.p_rrd_off_max*(1-pr.u_on) for pr in sc_data.prarray; lcon = fill(-Inf, size(sc_data.prarray)))
    c105_cs = constraint(core, p_jt_rrd_off_cs[cs.j_cs, cs.t] - cs.p_rrd_off_max*(1-cs.u_on) for cs in sc_data.csarray; lcon = fill(-Inf, size(sc_data.csarray)))
    #These constraints could be removed and the variables removed to simplify other constraints. However, they are kept for continuity
    c106 = constraint(core, p_jt_rrd_off_pr[pr.j_pr, pr.t] for pr in sc_data.prarray)
    c107 = constraint(core, p_jt_nsc_cs[cs.j_cs, cs.t] for cs in sc_data.csarray)
    c108 = constraint(core, p_jt_rru_off_cs[cs.j_cs, cs.t] for cs in sc_data.csarray)

    #4.6.7 Relative reserve limits, based on headroom to max/min, producing devices
    c109 = constraint(core, p_jt_on_pr[pr.j_pr, pr.t] + p_jt_rgu_pr[pr.j_pr, pr.t] + p_jt_scr_pr[pr.j_pr, pr.t] + p_jt_rru_on_pr[pr.j_pr, pr.t] - pr.p_max*pr.u_on for pr in sc_data.prarray; lcon = fill(-Inf, size(sc_data.prarray)))
    c110 = constraint(core, p_jt_on_pr[pr.j_pr, pr.t] - p_jt_rgd_pr[pr.j_pr, pr.t] - p_jt_rrd_on_pr[pr.j_pr, pr.t] - pr.p_min*pr.u_on for pr in sc_data.prarray; ucon = fill(Inf, size(sc_data.prarray)))
    c111 = constraint(core, p_jt_su_pr[pr.j_pr, pr.t] + p_jt_sd_pr[pr.j_pr, pr.t] + p_jt_nsc_pr[pr.j_pr, pr.t] + p_jt_rru_off_pr[pr.j_pr, pr.t] - pr.p_max*(1-pr.u_on) for pr in sc_data.prarray; lcon = fill(-Inf, size(sc_data.prarray)))
    c112 = constraint(core, q_jt_pr[pr.j_pr, pr.t] + q_jt_qru_pr[pr.j_pr, pr.t] - pr.q_max*(pr.u_on + pr.sum2_T_supc_pr_jt + pr.sum2_T_supc_pr_jt) for pr in sc_data.prarray; lcon = fill(-Inf, size(sc_data.prarray)))
    c113 = constraint(core, q_jt_pr[pr.j_pr, pr.t] - q_jt_qrd_pr[pr.j_pr, pr.t] - pr.q_min*(pr.u_on + pr.sum2_T_supc_pr_jt + pr.sum2_T_supc_pr_jt) for pr in sc_data.prarray; ucon = fill(Inf, size(sc_data.prarray)))

    c114 = constraint(core, q_jt_pr[pr.j_pr, pr.t] + q_jt_qru_pr[pr.j_pr, pr.t] - pr.q_max_p0*(pr.u_on + pr.sum2_T_supc_pr_jt + pr.sum2_T_sdpc_pr_jt) - pr.beta_max*p_jt_pr[pr.j_pr, pr.t] for pr in sc_data.prarray_pqbounds;
    lcon = fill(-Inf, size(sc_data.prarray_pqbounds)))
    c115 = constraint(core, q_jt_pr[pr.j_pr, pr.t] - q_jt_qrd_pr[pr.j_pr, pr.t] - pr.q_min_p0*(pr.u_on + pr.sum2_T_supc_pr_jt + pr.sum2_T_sdpc_pr_jt) - pr.beta_min*p_jt_pr[pr.j_pr, pr.t] for pr in sc_data.prarray_pqbounds;
    ucon = fill(Inf, size(sc_data.prarray_pqbounds)))
    c116 = constraint(core, q_jt_pr[pr.j_pr, pr.t] - pr.q_p0*(pr.u_on + pr.sum2_T_supc_pr_jt + pr.sum2_T_sdpc_pr_jt) - pr.beta*p_jt_pr[pr.j_pr, pr.t] for pr in sc_data.prarray_pqe)
    c117 = constraint(core, q_jt_qru_pr[pr.j_pr, pr.t] for pr in sc_data.prarray_pqe)
    c118 = constraint(core, q_jt_qrd_pr[pr.j_pr, pr.t] for pr in sc_data.prarray_pqe)

    #4.6.8 Relative reserve limits, based on headroom to max/min, consuming devices
    c119 = constraint(core, p_jt_on_cs[cs.j_cs, cs.t] + p_jt_rgd_cs[cs.j_cs, cs.t] + p_jt_rrd_on_cs[cs.j_cs, cs.t] - cs.p_max*cs.u_on for cs in sc_data.csarray; lcon = fill(-Inf, size(sc_data.csarray)))
    c120 = constraint(core, p_jt_on_cs[cs.j_cs, cs.t] - p_jt_rgu_cs[cs.j_cs, cs.t] - p_jt_scr_cs[cs.j_cs, cs.t] - p_jt_rru_on_cs[cs.j_cs, cs.t] -cs.p_min*cs.u_on for cs in sc_data.csarray; ucon = fill(Inf, size(sc_data.csarray)))
    c121 = constraint(core, p_jt_su_cs[cs.j_cs, cs.t] + p_jt_sd_cs[cs.j_cs, cs.t] + p_jt_rrd_off_cs[cs.j_cs, cs.t] - cs.p_max*(1-cs.u_on) for cs in sc_data.csarray; lcon = fill(-Inf, size(sc_data.csarray)))
    c122 = constraint(core, q_jt_cs[cs.j_cs, cs.t] + q_jt_qrd_cs[cs.j_cs, cs.t] - cs.q_max*(cs.u_on + cs.sum2_T_supc_cs_jt + cs.sum2_T_sdpc_cs_jt) for cs in sc_data.csarray; lcon = fill(-Inf, size(sc_data.csarray)))
    c123 = constraint(core, q_jt_cs[cs.j_cs, cs.t] - q_jt_qru_cs[cs.j_cs, cs.t] - cs.q_min*(cs.u_on + cs.sum2_T_supc_cs_jt + cs.sum2_T_sdpc_cs_jt) for cs in sc_data.csarray; ucon = fill(Inf, size(sc_data.csarray)))
    c124 = constraint(core, q_jt_cs[cs.j_cs, cs.t] + q_jt_qrd_cs[cs.j_cs, cs.t] - cs.q_max_p0*(cs.u_on + cs.sum2_T_supc_cs_jt + cs.sum2_T_sdpc_cs_jt) - cs.beta_max*p_jt_cs[cs.j_cs, cs.t] for cs in sc_data.csarray_pqbounds; lcon = fill(-Inf, size(sc_data.csarray_pqbounds)))
    c125 = constraint(core, q_jt_cs[cs.j_cs, cs.t] - q_jt_qru_cs[cs.j_cs, cs.t] - cs.q_min_p0*(cs.u_on + cs.sum2_T_supc_cs_jt + cs.sum2_T_sdpc_cs_jt) - cs.beta_min*p_jt_cs[cs.j_cs, cs.t] for cs in sc_data.csarray_pqbounds; ucon = fill(Inf, size(sc_data.csarray_pqbounds)))
    c126 = constraint(core, q_jt_cs[cs.j_cs, cs.t] - cs.q_p0*(cs.u_on + cs.sum2_T_supc_cs_jt + cs.sum2_T_sdpc_cs_jt) - cs.beta*p_jt_cs[cs.j_cs, cs.t] for cs in sc_data.csarray_pqe)
    c127 = constraint(core, q_jt_qru_cs[cs.j_cs, cs.t] for cs in sc_data.csarray_pqe)
    c128 = constraint(core, q_jt_qrd_cs[cs.j_cs, cs.t] for cs in sc_data.csarray_pqe)

    #4.6.9 Energy cost and value
    c130_pr = constraint(core, p_jt_pr[pr.j_pr, pr.t] for pr in sc_data.prarray)
    c130_pr_a = constraint!(core, c130_pr, pr.j_pr + L_J_pr*(pr.t-1) => -p_jtm_pr[pr.flat_k] for pr in sc_data.p_jtm_flattened_pr)
    c130_cs = constraint(core, p_jt_cs[cs.j_cs, cs.t] for cs in sc_data.csarray)
    c130_cs_a = constraint!(core, c130_cs, cs.j_cs + L_J_cs*(cs.t-1) => -p_jtm_cs[cs.flat_k] for cs in sc_data.p_jtm_flattened_cs)
    c131_pr = constraint(core, z_jt_en_pr[pr.j_pr, pr.t]/pr.dt for pr in sc_data.prarray)
    c131_pr_a = constraint!(core, c131_pr, pr.j_pr + L_J_pr*(pr.t-1) => -pr.c_en*p_jtm_pr[pr.flat_k] for pr in sc_data.p_jtm_flattened_pr)
    c131_cs = constraint(core, z_jt_en_cs[cs.j_cs, cs.t]/cs.dt for cs in sc_data.csarray)
    c131_cs_a = constraint!(core, c131_cs, cs.j_cs + L_J_cs*(cs.t-1) => -cs.c_en*p_jtm_cs[cs.flat_k] for cs in sc_data.p_jtm_flattened_cs)

    #4.7 Shunt devices
    c132 = constraint(core, p_jt_sh[sh.j_sh, sh.t] - g_jt_sh[sh.j_sh, sh.t]*v_it[sh.bus, sh.t]^2 for sh in sc_data.shuntarray)
    c133 = constraint(core, q_jt_sh[sh.j_sh, sh.t] + b_jt_sh[sh.j_sh, sh.t]*v_it[sh.bus, sh.t]^2 for sh in sc_data.shuntarray)
    c134 = constraint(core, g_jt_sh[sh.j_sh, sh.t] - sh.g_sh*sh.u_sh for sh in sc_data.shuntarray)
    c135 = constraint(core, b_jt_sh[sh.j_sh, sh.t] - sh.b_sh*sh.u_sh for sh in sc_data.shuntarray)
    #Assume (136-137) properly handled in uc solution

    #4.8.1 AC branch flow limits and penalties
    #AC branches split into ln and xf
    c139_ln = constraint(core, z_jt_s_ln[ln.j_ln, ln.t] - ln.dt*c_s*s_jt_plus_ln[ln.j_ln, ln.t] for ln in sc_data.aclbrancharray)
    c139_xf = constraint(core, z_jt_s_xf[xf.j_xf, xf.t] -xf.dt*c_s*s_jt_plus_xf[xf.j_xf, xf.t] for xf in sc_data.acxbrancharray)
    c140_ln = constraint(core, (p_jt_fr_ln[ln.j_ln, ln.t]^2 + q_jt_fr_ln[ln.j_ln, ln.t]^2)^.5 - ln.s_max - s_jt_plus_ln[ln.j_ln, ln.t] for ln in sc_data.aclbrancharray;
                        lcon = fill(-Inf, size(sc_data.aclbrancharray)))
    c140_xf = constraint(core, (p_jt_fr_xf[xf.j_xf, xf.t]^2 + q_jt_fr_xf[xf.j_xf, xf.t]^2)^.5 - xf.s_max - s_jt_plus_xf[xf.j_xf, xf.t] for xf in sc_data.acxbrancharray;
                        lcon = fill(-Inf, size(sc_data.acxbrancharray)))
    c141_ln = constraint(core, (p_jt_to_ln[ln.j_ln, ln.t]^2 + q_jt_to_ln[ln.j_ln, ln.t]^2)^.5 - ln.s_max - s_jt_plus_ln[ln.j_ln, ln.t] for ln in sc_data.aclbrancharray;
                        lcon = fill(-Inf, size(sc_data.aclbrancharray)))
    c141_xf = constraint(core, (p_jt_to_xf[xf.j_xf, xf.t]^2 + q_jt_to_xf[xf.j_xf, xf.t]^2)^.5 - xf.s_max - s_jt_plus_xf[xf.j_xf, xf.t] for xf in sc_data.acxbrancharray;
                        lcon = fill(-Inf, size(sc_data.acxbrancharray)))

    #4.8.2 AC branch controls
    c142 = constraint(core, φ_jt_ln[ln.j_ln, ln.t] for ln in sc_data.aclbrancharray)
    c143 = constraint(core, τ_jt_ln[ln.j_ln, ln.t] - 1 for ln in sc_data.aclbrancharray)
    c144 = constraint(core, φ_jt_xf[xf.j_xf, xf.t] - xf.phi_o for xf in sc_data.fpdarray)
    c145 = constraint(core, τ_jt_xf[xf.j_xf, xf.t] - xf.tau_o for xf in sc_data.fwrarray)
    c146 = constraint(core, φ_jt_xf[xf.j_xf, xf.t] for xf in sc_data.vpdarray;
                lcon = [xf.phi_min for xf in sc_data.vpdarray], ucon = [xf.phi_max for xf in sc_data.vpdarray])
    c147 = constraint(core, τ_jt_xf[xf.j_xf, xf.t] for xf in sc_data.vwrarray;
                lcon = [xf.tau_min for xf in sc_data.vwrarray], ucon = [xf.tau_max for xf in sc_data.vwrarray])

    #4.8.3 AC branch flows
    
    c148_ln = constraint(core, -p_jt_fr_ln[ln.j_ln, ln.t] + ln.u_on*((ln.g_sr+ln.g_fr)*v_it[ln.fr_bus, ln.t]^2/(τ_jt_ln[ln.j_ln, ln.t]^2) + (-ln.g_sr*cos(θ_it[ln.fr_bus, ln.t] - θ_it[ln.to_bus, ln.t] - φ_jt_ln[ln.j_ln, ln.t])
     - ln.b_sr*sin(θ_it[ln.fr_bus, ln.t] - θ_it[ln.to_bus, ln.t] - φ_jt_ln[ln.j_ln, ln.t]))*v_it[ln.fr_bus, ln.t]*v_it[ln.to_bus, ln.t]/τ_jt_ln[ln.j_ln, ln.t]) for ln in sc_data.aclbrancharray)
    c148_xf = constraint(core, -p_jt_fr_xf[xf.j_xf, xf.t] + xf.u_on*((xf.g_sr+xf.g_fr)*v_it[xf.fr_bus, xf.t]^2/(τ_jt_xf[xf.j_xf, xf.t]^2) + (-xf.g_sr*cos(θ_it[xf.fr_bus, xf.t] - θ_it[xf.to_bus, xf.t] - φ_jt_xf[xf.j_xf, xf.t])
     - xf.b_sr*sin(θ_it[xf.fr_bus, xf.t] - θ_it[xf.to_bus, xf.t] - φ_jt_xf[xf.j_xf, xf.t]))*v_it[xf.fr_bus, xf.t]*v_it[xf.to_bus, xf.t]/τ_jt_xf[xf.j_xf, xf.t]) for xf in sc_data.acxbrancharray)
    c149_ln = constraint(core, -q_jt_fr_ln[ln.j_ln, ln.t] + ln.u_on*((-ln.b_sr - ln.b_fr - ln.b_ch/2)*v_it[ln.fr_bus, ln.t]^2/(τ_jt_ln[ln.j_ln, ln.t]^2) + (ln.b_sr*cos(θ_it[ln.fr_bus, ln.t] - θ_it[ln.to_bus, ln.t] - φ_jt_ln[ln.j_ln, ln.t]) 
    - ln.g_sr*sin(θ_it[ln.fr_bus, ln.t] - θ_it[ln.to_bus, ln.t] - φ_jt_ln[ln.j_ln, ln.t]))*v_it[ln.fr_bus, ln.t]*v_it[ln.to_bus, ln.t]/τ_jt_ln[ln.j_ln, ln.t]) for ln in sc_data.aclbrancharray)
    c149_xf = constraint(core, -q_jt_fr_xf[xf.j_xf, xf.t] + xf.u_on*((-xf.b_sr - xf.b_fr - xf.b_ch/2)*v_it[xf.fr_bus, xf.t]^2/(τ_jt_xf[xf.j_xf, xf.t]^2) + (xf.b_sr*cos(θ_it[xf.fr_bus, xf.t] - θ_it[xf.to_bus, xf.t] - φ_jt_xf[xf.j_xf, xf.t]) 
    - xf.g_sr*sin(θ_it[xf.fr_bus, xf.t] - θ_it[xf.to_bus, xf.t] - φ_jt_xf[xf.j_xf, xf.t]))*v_it[xf.fr_bus, xf.t]*v_it[xf.to_bus, xf.t]/τ_jt_xf[xf.j_xf, xf.t]) for xf in sc_data.acxbrancharray)
    c150_ln = constraint(core, -p_jt_to_ln[ln.j_ln, ln.t] + ln.u_on*((ln.g_sr+ln.g_to)*v_it[ln.to_bus, ln.t]^2 + (-ln.g_sr*cos(θ_it[ln.fr_bus, ln.t] - θ_it[ln.to_bus, ln.t] - φ_jt_ln[ln.j_ln, ln.t]) 
    + ln.b_sr*sin(θ_it[ln.fr_bus, ln.t] - θ_it[ln.to_bus, ln.t] - φ_jt_ln[ln.j_ln, ln.t]))*v_it[ln.fr_bus, ln.t]*v_it[ln.to_bus, ln.t]/τ_jt_ln[ln.j_ln, ln.t]) for ln in sc_data.aclbrancharray)
    c150_xf = constraint(core, -p_jt_to_xf[xf.j_xf, xf.t] + xf.u_on*((xf.g_sr+xf.g_to)*v_it[xf.to_bus, xf.t]^2 + (-xf.g_sr*cos(θ_it[xf.fr_bus, xf.t] - θ_it[xf.to_bus, xf.t] - φ_jt_xf[xf.j_xf, xf.t]) 
    + xf.b_sr*sin(θ_it[xf.fr_bus, xf.t] - θ_it[xf.to_bus, xf.t] - φ_jt_xf[xf.j_xf, xf.t]))*v_it[xf.fr_bus, xf.t]*v_it[xf.to_bus, xf.t]/τ_jt_xf[xf.j_xf, xf.t]) for xf in sc_data.acxbrancharray)
    c151_ln = constraint(core, -q_jt_to_ln[ln.j_ln, ln.t] + ln.u_on*((-ln.b_sr-ln.b_to-ln.b_ch/2)*v_it[ln.to_bus, ln.t]^2 + (ln.b_sr*cos(θ_it[ln.fr_bus, ln.t] - θ_it[ln.to_bus, ln.t] - φ_jt_ln[ln.j_ln, ln.t]) 
    + ln.g_sr*sin(θ_it[ln.fr_bus, ln.t] - θ_it[ln.to_bus, ln.t] - φ_jt_ln[ln.j_ln, ln.t]))*v_it[ln.fr_bus, ln.t]*v_it[ln.to_bus, ln.t]/τ_jt_ln[ln.j_ln, ln.t]) for ln in sc_data.aclbrancharray)
    c151_xf = constraint(core, -q_jt_to_xf[xf.j_xf, xf.t] + xf.u_on*((-xf.b_sr-xf.b_to-xf.b_ch/2)*v_it[xf.to_bus, xf.t]^2 + (xf.b_sr*cos(θ_it[xf.fr_bus, xf.t] - θ_it[xf.to_bus, xf.t] - φ_jt_xf[xf.j_xf, xf.t]) 
    + xf.g_sr*sin(θ_it[xf.fr_bus, xf.t] - θ_it[xf.to_bus, xf.t] - φ_jt_xf[xf.j_xf, xf.t]))*v_it[xf.fr_bus, xf.t]*v_it[xf.to_bus, xf.t]/τ_jt_xf[xf.j_xf, xf.t]) for xf in sc_data.acxbrancharray)
    

    #4.8.4 DC lines
    c156 = constraint(core, p_jt_fr_dc[dc.j_dc, dc.t] + p_jt_to_dc[dc.j_dc, dc.t] for dc in sc_data.dclinearray)


    model = ExaModel(core; kwargs...)

    obj_vars = (
        z_w_en_max_pr = z_w_en_max_pr,
        z_w_en_max_cs = z_w_en_max_cs,
        z_w_en_min_pr = z_w_en_min_pr,
        z_w_en_min_cs = z_w_en_min_cs,
        
        #split z_jt_en and on into pr and cs
        z_jt_en_pr = z_jt_en_pr,
        z_jt_en_cs = z_jt_en_cs,

        z_it_p = z_it_p,
        z_it_q = z_it_q,

        #z_jt_s split into ln and xf
        z_jt_s_ln = z_jt_s_ln,
        z_jt_s_xf = z_jt_s_xf ,
        #z_jt rgu, rgd, scr, nsc, rru, rrd, qru, qrd split into pr and cs
        z_jt_rgu_pr = z_jt_rgu_pr, 
        z_jt_rgu_cs = z_jt_rgu_cs,
        z_jt_rgd_pr = z_jt_rgd_pr,
        z_jt_rgd_cs = z_jt_rgd_cs,
        z_jt_scr_pr = z_jt_scr_pr,
        z_jt_scr_cs = z_jt_scr_cs,
        z_jt_nsc_pr = z_jt_nsc_pr,
        z_jt_nsc_cs = z_jt_nsc_cs,
        z_jt_rru_pr = z_jt_rru_pr,
        z_jt_rru_cs = z_jt_rru_cs,
        z_jt_rrd_pr = z_jt_rrd_pr,
        z_jt_rrd_cs = z_jt_rrd_cs,
        z_jt_qru_pr = z_jt_qru_pr,
        z_jt_qru_cs = z_jt_qru_cs,
        z_jt_qrd_pr = z_jt_qrd_pr,
        z_jt_qrd_cs = z_jt_qrd_cs,

        z_nt_rgu = z_nt_rgu,
        z_nt_rgd = z_nt_rgd,
        z_nt_scr = z_nt_scr ,
        z_nt_nsc = z_nt_nsc,
        z_nt_rru = z_nt_rru ,
        z_nt_rrd = z_nt_rrd,
        z_nt_qru = z_nt_qru,
        z_nt_qrd = z_nt_qrd,
        )

    unincluded_obj = -(z_jt_on_pr + z_jt_on_cs + z_jt_su_pr + z_jt_su_cs + z_jt_su_ln + z_jt_su_xf
                    + z_jt_sd_pr + z_jt_sd_cs + z_jt_sd_ln + z_jt_sd_xf + z_jt_sus_pr + z_jt_sus_cs)
    return model, sc_data, obj_vars, unincluded_obj

end

