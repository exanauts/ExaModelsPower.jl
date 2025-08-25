module ExaModelsPower

import JLD2
import Downloads
import ExaModels: ExaCore, variable, constraint, ExaModel, objective, constraint!, convert_array, solution


include("parser.jl")
include("opf.jl")
include("scopf.jl")
include("mpopf.jl")
include("constraint.jl")
include("sc_parser.jl")

const NAMES = filter(names(@__MODULE__; all = true)) do x
    str = string(x)
    endswith(str, "model") && !startswith(str, "#")
end

for name in filter(names(@__MODULE__; all = true)) do x
    endswith(string(x), "model")
end
    @eval export $name
end
    
function __init__()
    if haskey(ENV, "EXA_MODELS_DEPOT")
        global TMPDIR = ENV["EXA_MODELS_DEPOT"]
    else
        global TMPDIR = joinpath(@__DIR__,"..","data")
        mkpath(TMPDIR)
    end
end

end # module ExaModelsExamples
