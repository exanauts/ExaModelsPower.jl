using Pkg
Pkg.activate("..")  # assumes docs/ is a subfolder of the main project

using Documenter, ExaModelsPower, DocumenterCitations, Literate

if !(@isdefined _PAGES)
    const _PAGES = [
        "Introduction" => "Introduction.md",
        "Tutorial" => [
            "opf_demo.md",
            "mpopf_demo.md"
        ],
        "OPF Formulations" => "opfs_doc.md",
        "API Manual" => "core.md"
    ]
end

if !(@isdefined _JL_FILENAMES)
    const _JL_FILENAMES = [
        "opf_demo.jl",
        "mpopf_demo.jl"
    ]
end

for jl_filename in _JL_FILENAMES

    Literate.markdown(
        joinpath(@__DIR__, "src", jl_filename),
        joinpath(@__DIR__, "src");
        documenter = true,
        execute = true,
    )

end


makedocs(;
    sitename = "My Documentation",
    modules = [ExaModelsPower],
    remotes = nothing,
    pages = _PAGES
)