using ExaModelsPower, MadNLP, ExaModels, MadNLPHSL, NLPModelsIpopt, NLPModels

model_cpu, ~ = opf_model("pglib_opf_case3_lmbd";)
~ = ipopt(model_cpu, tol = 1e-4, max_iter = 3, dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = 1e-4, linear_solver = "ma86")

m_cpu, v_cpu, c_cpu = opf_model("pglib_opf_case78484_epigrids")
result_ma86 = ipopt(m_cpu, tol = 1e-4, max_wall_time = Float64(900), dual_inf_tol=Float64(10000), constr_viol_tol=Float64(10000), compl_inf_tol=Float64(10000), bound_relax_factor = 1e-4, linear_solver = "ma86", honor_original_bounds = "no", print_timing_statistics = "yes", output_file = "ipopt_output")
