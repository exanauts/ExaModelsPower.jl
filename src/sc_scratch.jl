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
    L_J_cspr = length(data.sdd_ts_lookup)
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
                (j=j, j_dc = j_dc, pdc_max=pdc_max, qdc_fr_min=qdc_fr_min, qdc_to_min=qdc_to_min, qdc_fr_max=qdc_fr_max, qdc_to_max=qdc_to_max, to_bus=to_bus, fr_bus=fr_bus)
            end for val in values(data.dc_line_lookup)

        ]
    


    )
    return sc_data
end