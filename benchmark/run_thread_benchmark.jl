using ExaModelsPower, MadNLP, ExaModels, MadNLPHSL, NLPModelsIpopt, NLPModels, Plots, LinearAlgebra

BLAS.set_num_threads(1);
model_cpu, ~ = opf_model("pglib_opf_case3_lmbd";)
~ = madnlp(model_cpu, tol = 1e-4,
            kkt_system=MadNLP.SparseCondensedKKTSystem, equality_treatment=MadNLP.RelaxEquality, 
            fixed_variable_treatment=MadNLP.RelaxBound, dual_initialized=true,
            linear_solver=MadNLPHSL.Ma86Solver, ma86_num_threads=2)

m_cpu, v_cpu, c_cpu = opf_model("pglib_opf_case78484_epigrids")

threads = [12, 16, 20, 28]
solve_times = Float64[]

for thread in threads
    result_ma86 = madnlp(m_cpu, tol = 1e-4,
            kkt_system=MadNLP.SparseCondensedKKTSystem, equality_treatment=MadNLP.RelaxEquality, 
            fixed_variable_treatment=MadNLP.RelaxBound, dual_initialized=true,
            linear_solver=MadNLPHSL.Ma86Solver, ma86_num_threads=thread)
    push!(solve_times, result_ma86.counters.total_time)
end

# Plot Thread Count vs Solve Time
plot(
    threads, solve_times,
    marker = :o,
    xlabel = "Number of Threads",
    ylabel = "Solve Time (s)",
    title = "MA86 Threads vs Solve Time",
    legend = false,
)

savefig("thread_vs_solve_time_madnlp.png")
