#this code only needs to execute once to output solution file

#=Pkg.add(url = "https://github.com/lanl-ansi/GOC3Benchmark.jl.git")
path = pathof(GOC3Benchmark)  # gives the path to GOC3Benchmark.jl
include(joinpath(dirname(path), "..", "MyJulia1.jl"))
MyJulia1("data/C3E4N00073D1_scenario_303.json", 600, 1, "C3E4N00073", 1)=#




using GOC3Benchmark, JSON

#this solution found using the ARPA-E benchmark code
uc_data = JSON.parsefile("solution.json")
data = GOC3Benchmark.get_data_from_file("data/C3E4N00073D1_scenario_303.json")

function sc_parser(data)
    L_J_xf = length(data.twt_lookup)
    L_J_ln = length(data.ac_line_lookup)
    L_J_ac = L_J_ln + L_J_xf
    L_J_dc = length(data.dc_line_lookup)
    L_J_br =  L_J_ac + L_J_dc
    L_J_cs = length(data.sdd_ids_consumer)
    L_J_pr = length(data.sdd_ids_producer)
    L_J_cspr = L_J_cs + L_J_pr
    L_J_sh = length(data.shunt_lookup)
    
    sc_data = (
        bus = [
            begin
                i = parse(Int, match(r"\d+", val["uid"]).match)+1
                uid = val["uid"]
                v_min = val["vm_lb"]
                v_max = val["vm_ub"]
                (i = i, uid = uid, v_min = v_min, v_max = v_max)
            end for val in values(data.bus_lookup)
        ],

        shunt = [
            begin
                j = parse(Int, match(r"\d+", val["uid"]).match) + L_J_br + L_J_cspr+1
                j_sh = parse(Int, match(r"\d+", val["uid"]).match)+1
                uid = val["uid"]
                bus = parse(Int, match(r"\d+", val["bus"]).match)+1
                g_sh = val["gs"]
                b_sh = val["bs"]
                (j = j, uid = uid, g_sh = g_sh, b_sh = b_sh)
            end for val in values(data.shunt_lookup)
        ],

        ac_branch = vcat(
            # AC lines
            [
                begin
                    j = parse(Int, match(r"\d+", val["uid"]).match)+1
                    j_ac = j
                    j_ln = j
                    uid = val["uid"]
                    to_bus = parse(Int, match(r"\d+", val["to_bus"]).match)+1
                    fr_bus = parse(Int, match(r"\d+", val["fr_bus"]).match)+1
                    c_su = val["connection_cost"]
                    c_sd = val["disconnection_cost"]
                    s_max = val["mva_ub_nom"]
                    r = val["r"]
                    x = val["x"]
                    g_sr = r / (x^2 + r^2)
                    b_sr = -x / (x^2 + r^2)
                    b_ch = val["b"]
                    if val["additional_shunt"] == 1
                        g_fr = val["g_fr"]
                        g_to = val["g_to"]
                        b_fr = val["b_fr"]
                        b_to = val["b_to"]
                    else
                        g_fr = 0
                        g_to = 0
                        b_fr = 0
                        b_to = 0
                    end
                    (j = j, j_ac = j_ac, j_ln = j_ln, uid = uid, to_bus = to_bus, fr_bus = fr_bus, c_su = c_su, c_sd = c_sd, s_max = s_max, 
                    g_sr = g_sr, b_sr = b_sr, b_ch = b_ch, g_fr = g_fr, g_to = g_to, b_fr = b_fr, b_to = b_to)
                end for val in values(data.ac_line_lookup)
            ],
            
            # Transformers
            [
                begin
                    j_xf = parse(Int, match(r"\d+", val["uid"]).match)+1
                    j = j_xf + L_J_ln
                    j_ac = j
                    uid = val["uid"]
                    to_bus = parse(Int, match(r"\d+", val["to_bus"]).match)+1
                    fr_bus = parse(Int, match(r"\d+", val["fr_bus"]).match)+1
                    c_su = val["connection_cost"]
                    c_sd = val["disconnection_cost"]
                    s_max = val["mva_ub_nom"]
                    r = val["r"]
                    x = val["x"]
                    g_sr = r / (x^2 + r^2)
                    b_sr = -x / (x^2 + r^2)
                    b_ch = val["b"]
                    if val["additional_shunt"] == 1
                        g_fr = val["g_fr"]
                        g_to = val["g_to"]
                        b_fr = val["b_fr"]
                        b_to = val["b_to"]
                    else
                        g_fr = 0
                        g_to = 0
                        b_fr = 0
                        b_to = 0
                    end
                    (j = j, j_ac = j_ac, j_xf=j_xf, uid = uid, to_bus = to_bus, fr_bus = fr_bus, c_su = c_su, c_sd = c_sd, s_max = s_max, 
                    g_sr = g_sr, b_sr = b_sr, b_ch = b_ch, g_fr = g_fr, g_to = g_to, b_fr = b_fr, b_to = b_to)                
                end for val in values(data.twt_lookup)
            ]
        ),
        #Variable phase difference
        vpd = [
            begin
                j = parse(Int, match(r"\d+", val["uid"]).match) + L_J_ln+1
                j_ac = j
                phi_min = val["ta_lb"]
                phi_max = val["ta_ub"]
                (j = j, j_ac = j_ac, phi_min = phi_min, phi_max = phi_max)
            end for val in values(data.twt_lookup) if val["ta_lb"] < val["ta_ub"]
        ], 
        #Fixed phase difference
        fpd = [
            begin
                j = parse(Int, match(r"\d+", val["uid"]).match) + L_J_ln+1
                j_ac = j
                phi_o = val["initial_status"]["ta"]
                (j = j, j_ac = j_ac, phi_o = phi_o)
            end for val in values(data.twt_lookup) if val["ta_lb"] >= val["ta_ub"]
        ],
        #Variable winding ratio
        vwd = [
            begin
                j = parse(Int, match(r"\d+", val["uid"]).match) + L_J_ln+1
                j_ac = j
                tau_min = val["tm_lb"]
                tau_max = val["tm_ub"]
                (j=j, j_ac=j_ac, tau_min=tau_min, tau_max=tau_max)
            end for val in values(data.twt_lookup) if val["tm_lb"] < val["tm_ub"]
        ],
        #Fixed winding ratio
        fwd = [
            begin
                j = parse(Int, match(r"\d+", val["uid"]).match) + L_J_ln+1
                j_ac = j
                tau_o = val["initial_status"]["tm"]
                (j=j, j_ac=j_ac, tau_o=tau_o)
            end for val in values(data.twt_lookup) if val["tm_lb"] >= val["tm_ub"]

        ],

        dc_branch = [
            begin
                j = parse(Int, match(r"\d+", val["uid"]).match) + L_J_ac+1
                j_dc = parse(Int, match(r"\d+", val["uid"]).match)+1
                uid = val["uid"]
                pdc_max = val["pdc_ub"]
                qdc_fr_min = val["qdc_fr_lb"]
                qdc_to_min = val["qdc_to_lb"]
                qdc_fr_max = val["qdc_fr_ub"]
                qdc_to_max = val["qdc_to_ub"]
                to_bus = parse(Int, match(r"\d+", val["to_bus"]).match)
                fr_bus = parse(Int, match(r"\d+", val["fr_bus"]).match)
                (j=j, j_dc = j_dc, uid=uid, pdc_max=pdc_max, qdc_fr_min=qdc_fr_min, qdc_to_min=qdc_to_min, qdc_fr_max=qdc_fr_max, qdc_to_max=qdc_to_max, to_bus=to_bus, fr_bus=fr_bus)
            end for val in values(data.dc_line_lookup)

        ],

        prod_cons = vcat(
            # Producers
            [
                begin
                    ts_val = data.sdd_ts_lookup[key]
                    j = parse(Int, match(r"\d+", val["uid"]).match) + L_J_br + 1
                    j_prcs = parse(Int, match(r"\d+", val["uid"]).match) + 1
                    j_pr = parse(Int, match(r"\d+", val["uid"]).match) + 1
                    uid = val["uid"]
                    c_on = val["on_cost"]
                    c_su = val["startup_cost"]
                    c_sd = val["shutdown_cost"]
                    p_ru = val["p_ramp_up_ub"]
                    p_rd = val["p_ramp_down_ub"]
                    p_ru_su = val["p_startup_ramp_ub"]
                    p_rs_sd = val["p_shutdown_ramp_ub"]
                    c_rgu = ts_val["p_reg_res_up_cost"]
                    c_rgd = ts_val["p_reg_res_down_cost"]
                    c_scr = ts_val["p_syn_res_cost"]
                    c_nsc = ts_val["p_nsyn_res_cost"]
                    c_rru_on = ts_val["p_ramp_res_up_online_cost"]
                    c_rru_off = ts_val["p_ramp_res_up_offline_cost"]
                    c_rrd_on = ts_val["p_ramp_res_down_online_cost"]
                    c_rrd_off = ts_val["p_ramp_res_down_offline_cost"]
                    c_qru = ts_val["q_res_up_cost"]
                    c_qrd = ts_val["q_res_down_cost"]
                    p_rgu_max = val["p_reg_res_up_ub"]
                    p_rgd_max = val["p_reg_res_down_ub"]
                    p_scr_max = val["p_syn_res_ub"]
                    p_nsc_max = val["p_nsyn_res_ub"]
                    p_rru_on_max = val["p_ramp_res_up_online_ub"]
                    p_rru_off_max = val["p_ramp_res_up_offline_ub"]
                    p_rrd_on_max = val["p_ramp_res_down_online_ub"]
                    p_rrd_off_max = val["p_ramp_res_down_offline_ub"]
                    
                    (j = j, j_prcs = j_prcs, j_pr = j_pr, uid = uid, c_on = c_on, c_su = c_su, c_sd = c_sd, p_ru = p_ru, p_rd = p_rd, p_ru_su = p_ru_su, p_rs_sd = p_rs_sd, 
                    c_rgu = c_rgu, c_rg = c_rgd, c_scr = c_scr, c_nsc = c_nsc, c_rru_on = c_rru_on, c_rru_off = c_rru_off, c_rrd_on = c_rrd_on, c_rrd_off = c_rrd_off, 
                    c_qru = c_qru, c_qrd = c_qrd, p_rgu_max = p_rgu_max, p_rgd_max = p_rgd_max, p_scr_max = p_scr_max, p_nsc_max = p_nsc_max, p_rru_on_max = p_rru_on_max,
                    p_rru_off_max, p_rrd_on_max, p_rrd_off_max)
                end for (key, val) in data.sdd_lookup if val["device_type"] == "producer"
            ],
            
            #Consumers
            [
                begin
                    ts_val = data.sdd_ts_lookup[key]
                    j = parse(Int, match(r"\d+", val["uid"]).match) + L_J_br + 1
                    j_prcs = parse(Int, match(r"\d+", val["uid"]).match) + 1
                    j_cs = parse(Int, match(r"\d+", val["uid"]).match) + L_J_pr + 1
                    uid = val["uid"]
                    c_on = val["on_cost"]
                    c_su = val["startup_cost"]
                    c_sd = val["shutdown_cost"]
                    p_ru = val["p_ramp_up_ub"]
                    p_rd = val["p_ramp_down_ub"]
                    p_ru_su = val["p_startup_ramp_ub"]
                    p_rs_sd = val["p_shutdown_ramp_ub"]
                    c_rgu = ts_val["p_reg_res_up_cost"]
                    c_rgd = ts_val["p_reg_res_down_cost"]
                    c_scr = ts_val["p_syn_res_cost"]
                    c_nsc = ts_val["p_nsyn_res_cost"]
                    c_rru_on = ts_val["p_ramp_res_up_online_cost"]
                    c_rru_off = ts_val["p_ramp_res_up_offline_cost"]
                    c_rrd_on = ts_val["p_ramp_res_down_online_cost"]
                    c_rrd_off = ts_val["p_ramp_res_down_offline_cost"]
                    c_qru = ts_val["q_res_up_cost"]
                    c_qrd = ts_val["q_res_down_cost"]
                    p_rgu_max = val["p_reg_res_up_ub"]
                    p_rgd_max = val["p_reg_res_down_ub"]
                    p_scr_max = val["p_syn_res_ub"]
                    p_nsc_max = val["p_nsyn_res_ub"]
                    p_rru_on_max = val["p_ramp_res_up_online_ub"]
                    p_rru_off_max = val["p_ramp_res_up_offline_ub"]
                    p_rrd_on_max = val["p_ramp_res_down_online_ub"]
                    p_rrd_off_max = val["p_ramp_res_down_offline_ub"]
                    
                    (j = j, j_prcs = j_prcs, j_cs = j_cs, uid = uid, c_on = c_on, c_su = c_su, c_sd = c_sd, p_ru = p_ru, p_rd = p_rd, p_ru_su = p_ru_su, p_rs_sd = p_rs_sd, 
                    c_rgu = c_rgu, c_rg = c_rgd, c_scr = c_scr, c_nsc = c_nsc, c_rru_on = c_rru_on, c_rru_off = c_rru_off, c_rrd_on = c_rrd_on, c_rrd_off = c_rrd_off, 
                    c_qru = c_qru, c_qrd = c_qrd, p_rgu_max = p_rgu_max, p_rgd_max = p_rgd_max, p_scr_max = p_scr_max, p_nsc_max = p_nsc_max, p_rru_on_max = p_rru_on_max,
                    p_rru_off_max, p_rrd_on_max, p_rrd_off_max)                
                end for (key, val) in data.sdd_lookup if val["device_type"] == "consumer"
            ]
        ),
        active_reserve = [
            begin
                ts_val = data.azr_ts_lookup[key]
                j = parse(Int, match(r"\d+", val["uid"]).match) + 1
                uid = val["uid"]
                c_rgu = val["REG_UP_vio_cost"]
                c_rgd = val["REG_DOWN_vio_cost"]
                c_scr = val["SYN_vio_cost"]
                c_nsc = val["NSYN_vio_cost"]
                c_rru = val["RAMPING_RESERVE_UP_vio_cost"]
                c_rrd = val["RAMPING_RESERVE_DOWN_vio_cost"]
                σ_rgu = val["REG_UP"]
                σ_rgd = val["REG_DOWN"]
                σ_scr = val["SYN"]
                σ_nsc = val["NSYN"]
                p_rru_min = ts_val["RAMPING_RESERVE_UP"]
                p_rrd_min = ts_val["RAMPING_RESERVE_DOWN"]
                (j=j, uid=uid, c_rgu=c_rgu, c_rgd=c_rgd, c_scr=c_scr, c_nsc=c_nsc, c_rru=c_rru, c_rrd=c_rrd, σ_rgu=σ_rgu, σ_rgd=σ_rgd, σ_scr=σ_scr, 
                σ_nsc=σ_nsc, p_rru_min=p_rru_min, p_rrd_min=p_rrd_min)
            end for (key, val) in data.azr_lookup

        ],
    

    )
    return sc_data
end


