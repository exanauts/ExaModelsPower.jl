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

    v_lvar = repeat([b.v_min for b in sc_data.bus], 1, L_T)
    v_uvar = repeat([b.v_max for b in sc_data.bus], 1, L_T)
    sc_data = convert_data(sc_data, backend)

    core = ExaCore(T; backend =backend)

    #variables are indexed j,t,k or j,t (t always second if present)

    b_jt_sh = variable(core, L_J_sh, L_T;)
    #Split e_w_plus into separate sets for W_en_min and W_en_max ad for pr, cs
    e_w_plus_min_pr = variable(core, L_W_en_min_pr;)
    e_w_plus_min_cs = variable(core, L_W_en_min_cs;)
    e_w_plus_max_pr = variable(core, L_W_en_max_pr;)
    e_w_plus_max_cs = variable(core, L_W_en_max_cs;)
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
    p_jtm = variable(core, length(sc_data.prcsflattened);)
    #to/from power split into ln, xf, and dc lines
    p_jt_fr_ln = variable(core, L_J_ln, L_T;)
    p_jt_fr_xf = variable(core, L_J_xf, L_T;)
    p_jt_fr_dc = variable(core, L_J_dc, L_T;)
    p_jt_to_ln = variable(core, L_J_ln, L_T;)
    p_jt_to_xf = variable(core, L_J_xf, L_T;)
    p_jt_to_dc = variable(core, L_J_dc, L_T;)
    q_jt_fr_ln = variable(core, L_J_ln, L_T;)
    q_jt_fr_xf = variable(core, L_J_xf, L_T;)
    q_jt_fr_dc = variable(core, L_J_dc, L_T;)
    q_jt_to_ln = variable(core, L_J_ln, L_T;)
    q_jt_to_xf = variable(core, L_J_xf, L_T;)
    q_jt_to_dc = variable(core, L_J_dc, L_T;)
    #p_jt rgu, rgd, scr, rru,on, rru,off, rrd,on, rrd,off and q_jt qru/qrd split into pr and cs
    p_jt_rgu_pr = variable(core, L_J_pr, L_T;)
    p_jt_rgu_cs = variable(core, L_J_cs, L_T;)
    p_jt_rgd_pr = variable(core, L_J_pr, L_T;)
    p_jt_rgd_cs = variable(core, L_J_cs, L_T;)
    p_jt_scr_pr = variable(core, L_J_pr, L_T;)
    p_jt_scr_cs = variable(core, L_J_cs, L_T;)
    p_jt_nsc_pr = variable(core, L_J_pr, L_T;)
    p_jt_nsc_cs = variable(core, L_J_cs, L_T;)
    p_jt_rru_on_pr = variable(core, L_J_pr, L_T;)
    p_jt_rru_on_cs = variable(core, L_J_cs, L_T;)
    p_jt_rru_off_pr = variable(core, L_J_pr, L_T;)
    p_jt_rru_off_cs = variable(core, L_J_cs, L_T;)
    p_jt_rrd_on_pr = variable(core, L_J_pr, L_T;)
    p_jt_rrd_on_cs = variable(core, L_J_cs, L_T;)
    p_jt_rrd_off_pr = variable(core, L_J_pr, L_T;)
    p_jt_rrd_off_cs = variable(core, L_J_cs, L_T;)
    q_jt_qru_pr = variable(core, L_J_pr, L_T;)
    q_jt_qru_cs = variable(core, L_J_cs, L_T;)
    q_jt_qrd_pr = variable(core, L_J_pr, L_T;)
    q_jt_qrd_cs = variable(core, L_J_cs, L_T;)

    
    p_nt_rgu_req = variable(core, Np, L_T;)
    p_nt_rgd_req = variable(core, Np, L_T;)
    p_nt_scr_req = variable(core, Np, L_T;)
    p_nt_nsc_req = variable(core, Np, L_T;)

    p_jt_pr_max = variable(core, L_T;)
    #Bounds from 4.3.1 Reserve shortfall domains (20-27)
    p_nt_rgu_plus = variable(core, Np, L_T; lvar = zeros(size(sc_data.preservearray)))
    p_nt_rgd_plus = variable(core, Np, L_T; lvar = zeros(size(sc_data.preservearray)))
    p_nt_scr_plus = variable(core, Np, L_T; lvar = zeros(size(sc_data.preservearray)))
    p_nt_nsc_plus = variable(core, Np, L_T; lvar = zeros(size(sc_data.preservearray)))
    p_nt_rru_plus = variable(core, Np, L_T; lvar = zeros(size(sc_data.preservearray)))
    p_nt_rrd_plus = variable(core, Np, L_T; lvar = zeros(size(sc_data.preservearray)))
    q_nt_qru_plus = variable(core, Nq, L_T; lvar = zeros(size(sc_data.qreservearray)))
    q_nt_qrd_plus = variable(core, Nq, L_T; lvar = zeros(size(sc_data.qreservearray)))

    p_t_sl = variable(core, L_T;)

    q_it = variable(core, I, L_T;)
    q_it_plus = variable(core, I, L_T;)

    #s_jt_plus split on ln and xf
    s_jt_plus_ln = variable(core, L_J_ln, L_T;)
    s_jt_lus_xf = variable(core, L_J_xf, L_T;)

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
    z_jt_on_pr = variable(core, L_J_pr, L_T;)
    z_jt_on_cs = variable(core, L_J_cs, L_T;)

    z_it_p = variable(core, I, L_T;)
    z_it_q = variable(core, I, L_T;)

    #z su/sd split into pr, cs, ln, xf
    z_jt_sd_pr = variable(core, L_J_pr, L_T;)
    z_jt_sd_cs = variable(core, L_J_cs, L_T;)
    z_jt_sd_ln = variable(core, L_J_ln, L_T;)
    z_jt_sd_xf = variable(core, L_J_xf, L_T;)
    z_jt_su_pr = variable(core, L_J_pr, L_T;)
    z_jt_su_cs = variable(core, L_J_cs, L_T;)
    z_jt_su_ln = variable(core, L_J_ln, L_T;)
    z_jt_su_xf = variable(core, L_J_xf, L_T;)
    #z_jt_sus split into pr and cs
    z_jt_sus_pr = variable(core, L_J_pr, L_T;)
    z_jt_sus_cs = variable(core, L_J_cs, L_T;)
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
    τ_jt_ln = variable(core, L_J_ln, L_T;)
    τ_jt_xf = variable(core, L_J_xf, L_T;)
    φ_jt_ln = variable(core, L_J_ln, L_T;)
    φ_jt_xf = variable(core, L_J_xf, L_T;)

    #skip overall obj contraints for now
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

    model = ExaModel(core; kwargs...)
    return model, sc_data
end

