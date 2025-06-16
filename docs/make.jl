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
        "API Manual" => "core.md",
        "References" => "ref.md"
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

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(;
    plugins = [bib],
    sitename = "ExaModelsPower.jl",
    modules = [ExaModelsPower],
    remotes = nothing,
    authors = "Sungho Shin, Sanjay Johnson",
    format = Documenter.HTML(
        assets = ["assets/citations.css"],
        prettyurls = true,
        sidebar_sitename = true,
        collapselevel = 1,
    ),
    pages = _PAGES,
    clean = false
)

deploydocs(repo = "github.com/exanauts/ExaModelsPower.jl.git"; push_preview = true)

#- The step-by-step tutorials of using ExaModelsPower.jl can be found in [OPF tutorial](@ref opf_demo) and [MPOPF tutorial](@ref mpopf_demo).