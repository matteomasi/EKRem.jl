# Documentation of SimEK, built with Documenter.jl. The example scripts are converted to
# documentation pages by Literate.jl, which runs them and writes their output and figures into
# the pages. A script is run again only when it, the data in its folder or the source code of
# SimEK have changed since the last build.
#
# Usage (from the repository root):
#     julia --project=docs docs/make.jl          # static site in docs/build
#     julia --project=docs docs/make.jl serve    # local server with live reload (LiveServer.jl)
#     julia --project=docs docs/make.jl rerun    # run all examples again

ENV["GKSwstype"] = "100"   # GR (Plots) without a display window

const ROOT = dirname(@__DIR__)
const SRC = joinpath(@__DIR__, "src")
const RERUN = "rerun" in ARGS

# LiveServer builds the documentation by running this script, serves the site and runs the
# script again whenever a page or a file in examples/ changes. The generated example pages
# are not watched. (The `literate_dir` option of servedocs is not used: it fails with Windows
# paths.)
if "serve" in ARGS
    using LiveServer
    filter!(!=("serve"), ARGS)
    servedocs(;
        foldername=@__DIR__,
        include_dirs=[joinpath(ROOT, "examples")],
        skip_dirs=[joinpath(SRC, "examples")],
        launch_browser=true,
    )
    exit()
end

using Documenter
using SimEK
using Literate
using Base64
using SHA

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
# code and the dependencies of SimEK
function fingerprint(script)
    data = joinpath(dirname(script), "data")
    files = [script; isdir(data) ? readdir(data; join=true) : String[]; joinpath(ROOT, "Project.toml")]
    append!(files, sort([joinpath(r, f) for (r, _, fs) in walkdir(joinpath(ROOT, "src")) for f in fs]))
    context = SHA.SHA1_CTX()
    foreach(f -> SHA.update!(context, read(f)), files)
    return bytes2hex(SHA.digest!(context))
end

# Only the last line of the progress report of each run is shown in the pages
SimEK.PROGRESS_INTERVAL[] = Inf

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

# Plots shows a figure as an HTML image with the PNG data encoded in base64. The figures are
# saved as PNG files next to the page and linked from it, which keeps the HTML pages small.
function figures_to_files(dir, name)
    pattern = r"```@raw html\n<img src=\"data:image/png;base64,([A-Za-z0-9+/=]+)\"[^>]*>\n```"
    return function (content)
        k = 0
        return replace(content, pattern => function (block)
            k += 1
            file = "$name-$k.png"
            write(joinpath(dir, file), base64decode(match(pattern, block)[1]))
            return "![]($file)"
        end)
    end
end

function example_page(script, title)
    path = joinpath(ROOT, "examples", script)
    dir = joinpath(SRC, "examples", dirname(script))
    name = splitext(basename(script))[1]
    page = joinpath(dir, name * ".md")
    stamp = page * ".sha1"
    key = fingerprint(path)
    if RERUN || !isfile(page) || !isfile(stamp) || read(stamp, String) != key
        @info "Running example $script"
        Literate.markdown(
            path,
            dir;
            documenter=true,
            execute=true,
            credit=false,
            preprocess=hide_regions,
            postprocess=figures_to_files(dir, name),
        )
        write(stamp, key)
    else
        @info "Example $script unchanged: not run again"
    end
    return title => joinpath("examples", dirname(script), name * ".md")
end
example_pages = [group => [example_page(s, t) for (s, t) in scripts] for (group, scripts) in example_groups]

makedocs(;
    modules=[SimEK],
    sitename="SimEK.jl",
    authors="Matteo Masi",
    repo=Remotes.GitHub("matteomasi", "SimEK.jl"),
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", "false") == "true", edit_link="main"),
    pages=[
        "Home" => "index.md",
        "Theory" => "theory.md",
        "User guide" => [
            "Installation" => "guide/installation.md",
            "Model set-up" => "guide/model.md",
            "Geochemical reactions" => "guide/reactions.md",
            "Simulation and results" => "guide/simulation.md",
        ],
        "Examples" => ["Overview" => "examples.md"; example_pages],
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
