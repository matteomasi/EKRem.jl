# Documentation of EKRem: Documenter.jl with the VitePress backend (DocumenterVitepress.jl).
# The example scripts are converted to documentation pages by Literate.jl, which runs them
# and writes their output and figures into the pages. A script is run again only when it,
# the data in its folder or the source code of EKRem have changed since the last build.
#
# Usage (from the repository root):
#     julia --project=docs docs/make.jl          # static site in docs/build/1
#     julia --project=docs docs/make.jl serve    # local server with live reload
#     julia --project=docs docs/make.jl rerun    # run all examples again
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
using SHA

const SERVE = "serve" in ARGS
const RERUN = "rerun" in ARGS
const ROOT = dirname(@__DIR__)
const SRC = joinpath(@__DIR__, "src")

# Example scripts by group; the paths are relative to examples/
example_groups = [
    "1D benchmarks" => [
        "benchmarks/1D/1D_trace_electromigration.jl" => "Electromigration and electroosmosis of a trace cation",
        "benchmarks/1D/1D_binary_fixed_potential.jl" => "Binary electrolyte, fixed voltage",
        "benchmarks/1D/1D_binary_fixed_current.jl" => "Binary electrolyte, fixed current",
        "benchmarks/1D/1D_benchmark_krabbenhoft.jl" => "Krabbenhøft and Krabbenhøft benchmark",
        "benchmarks/1D/1D_benchmark_phreeqc.jl" => "PHREEQC Ampere benchmark",
    ],
    "Other examples" => ["misc/1D_compartments_fixed_voltage.jl" => "1D setup with electrode compartments and controlled voltage"],
]

# Fingerprint of an example: the script, the files in the data/ folder next to it, the source
# code and the dependencies of EKRem
function fingerprint(script)
    data = joinpath(dirname(script), "data")
    files = [script; isdir(data) ? readdir(data; join=true) : String[]; joinpath(ROOT, "Project.toml")]
    append!(files, sort([joinpath(r, f) for (r, _, fs) in walkdir(joinpath(ROOT, "src")) for f in fs]))
    context = SHA.SHA1_CTX()
    foreach(f -> SHA.update!(context, read(f)), files)
    return bytes2hex(SHA.digest!(context))
end

# Only the last line of the progress report of each run is shown in the pages
EKRem.PROGRESS_INTERVAL[] = Inf

# The code between the lines `#hide-begin` and `#hide-end` of an example is run, but it is not
# shown in the page: its lines are marked with #hide, its comments become code comments, and
# its empty lines are removed. Lines marked with #src are left unchanged.
function hide_regions(content)
    lines = String[]
    hidden = false
    for line in split(replace(content, "\r\n" => "\n"), '\n')
        s = strip(line)
        if s == "#hide-begin" || s == "#hide-end"
            hidden = s == "#hide-begin"
        elseif !hidden || endswith(s, "#src")
            push!(lines, line)
        elseif !isempty(s)
            push!(lines, (startswith(s, "#") ? "#" : "") * line * " #hide")
        end
    end
    return join(lines, '\n')
end

function example_page(script, title)
    path = joinpath(ROOT, "examples", script)
    dir = joinpath(SRC, "examples", dirname(script))
    page = joinpath(dir, replace(basename(script), ".jl" => ".md"))
    stamp = page * ".sha1"
    key = fingerprint(path)
    if RERUN || !isfile(page) || !isfile(stamp) || read(stamp, String) != key
        @info "Running example $script"
        Literate.markdown(path, dir; documenter=true, execute=true, credit=false, preprocess=hide_regions)
        write(stamp, key)
    else
        @info "Example $script unchanged: not run again"
    end
    return title => joinpath("examples", replace(script, ".jl" => ".md"))
end
example_pages = [group => [example_page(s, t) for (s, t) in scripts] for (group, scripts) in example_groups]

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
