module ExaModelsPower

import JLD2
import Downloads
import ExaModels
import PowerModels

include("parser.jl")
include("opf.jl")
include("scopf.jl")
include("mpopf.jl")

export opf_model, scopf_model, mpopf_model
    
function __init__()
    if haskey(ENV, "EXA_MODELS_DEPOT")
        global TMPDIR = ENV["EXA_MODELS_DEPOT"]
    else
        global TMPDIR = joinpath(@__DIR__,"..","data")
        mkpath(TMPDIR)
    end
    PowerModels.silence()
end

end # module ExaModelsExamples
