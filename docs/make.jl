# Documentation of EKRem: Documenter.jl with the VitePress backend (DocumenterVitepress.jl).
# The example scripts are converted to documentation pages by Literate.jl and run during the
# build.
#
# Usage (from the repository root):
#     julia --project=docs docs/make.jl          # static site in docs/build/1
#     julia --project=docs docs/make.jl serve    # local server with live reload
#
# The VitePress packages (docs/package.json) are installed by npm during the first build,
# which requires internet access. Node.js is provided by DocumenterVitepress, except on
# Windows, where DocumenterVitepress does not run VitePress: there, the site is built or
# served by the npm scripts of docs/package.json, with a system installation of Node.js.

ENV["GKSwstype"] = "100"   # GR (Plots) without a display window

using Documenter
using DocumenterVitepress
using EKRem
using Literate

const SERVE = "serve" in ARGS
const ROOT = dirname(@__DIR__)
const SRC = joinpath(@__DIR__, "src")

# Example scripts → Markdown pages (lines ending with `#src` are omitted, lines starting with
# `#md` appear only in the pages)
examples = [
    "1D_trace_electromigration.jl" => "Electromigration and electroosmosis of a trace ion",
    "1D_binary_fixed_potential.jl" => "Binary electrolyte, fixed voltage",
    "1D_binary_fixed_current.jl" => "Binary electrolyte, fixed current",
    "1D_benchmark_krabbenhoft.jl" => "Krabbenhøft benchmark",
]
example_pages = map(examples) do (script, title)
    Literate.markdown(joinpath(ROOT, "examples", script), joinpath(SRC, "examples"); documenter=true, credit=false)
    return title => joinpath("examples", replace(script, ".jl" => ".md"))
end

# Only the last line of the progress report of each run is shown in the pages
EKRem.PROGRESS_INTERVAL[] = Inf

makedocs(;
    modules=[EKRem],
    sitename="EKRem.jl",
    authors="Matteo Masi",
    format=DocumenterVitepress.MarkdownVitepress(;
        repo="github.com/matteomasi/EKRem",
        devbranch="main",
        devurl="dev",
        description="Simulation of electrokinetic remediation in Julia",
        build_vitepress=!SERVE,
    ),
    pages=[
        "Home" => "index.md",
        "Theory" => "theory.md",
        "User guide" => [
            "Installation" => "guide/installation.md",
            "Model set-up" => "guide/model.md",
            "Geochemical reactions" => "guide/reactions.md",
            "Simulation and results" => "guide/simulation.md",
        ],
        "Examples" => example_pages,
        "API reference" => [
            "Models" => "reference/models.md",
            "Species and porous medium" => "reference/medium.md",
            "Potential and boundaries" => "reference/boundaries.md",
            "Geochemistry" => "reference/geochemistry.md",
        ],
        "Development" => ["Development set-up" => "development.md", "Private functions" => "private.md"],
    ],
    checkdocs=:exports,
)

function vitepress(command)
    return cd(() -> success(pipeline(`cmd /c npm run docs:$command`; stdout, stderr)), uppercasefirst(@__DIR__))
end

if Sys.iswindows()
    SERVE ? vitepress("dev") : (vitepress("build") || error("VitePress build failed"))
elseif SERVE
    DocumenterVitepress.dev_docs(joinpath(@__DIR__, "build"))
end
