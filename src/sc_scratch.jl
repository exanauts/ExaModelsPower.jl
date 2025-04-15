#this code only needs to execute once to output solution file

#=Pkg.add(url = "https://github.com/lanl-ansi/GOC3Benchmark.jl.git")
using GOC3Benchmark
path = pathof(GOC3Benchmark)  # gives the path to GOC3Benchmark.jl
include(joinpath(dirname(path), "..", "MyJulia1.jl"))
MyJulia1("data/C3E4N00073D1_scenario_303.json", 600, 1, "C3E4N00073", 1)=#


uc_data = JSON.parsefile("solution.json")