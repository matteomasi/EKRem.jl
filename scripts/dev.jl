# Development tasks of EKRem, called by the Makefile:
#
#     julia --project=. scripts/dev.jl <task> [argument]
#
# Tasks: help, install, examples, coverage, lint, format, format-check,
#        bump <patch|minor|major>, clean.
#
# The development tools (Aqua, Coverage, JuliaFormatter) are installed by `install` in the
# shared Julia environment `ekrem-dev`, which is added to the load path when a task needs
# them. They are therefore not dependencies of EKRem.

using Pkg

const ROOT = dirname(@__DIR__)
const DEV_ENV = "ekrem-dev"
const DEV_TOOLS = ["Aqua", "Coverage", "JuliaFormatter"]
const CODE_DIRS = ["src", "test", "scripts", "examples"]

const TASK = isempty(ARGS) ? "help" : ARGS[1]
const TASK_ARG = length(ARGS) >= 2 ? ARGS[2] : ""

# Packages are loaded at top level, so that they are visible to the functions called below.
if TASK in ("lint", "coverage", "format", "format-check")
    push!(LOAD_PATH, "@" * DEV_ENV)
    try
        TASK == "lint" && @eval using Aqua, EKRem
        TASK == "coverage" && @eval using Coverage
        TASK in ("format", "format-check") && @eval using JuliaFormatter
    catch err
        err isa ArgumentError || rethrow()
        println(stderr, "Development tools not found. Run `make install-dev` first.")
        exit(1)
    end
end

code_dirs() = [joinpath(ROOT, d) for d in CODE_DIRS if isdir(joinpath(ROOT, d))]

"Print the Makefile targets that carry a `## description` comment."
function help()
    println("EKRem Development Commands:\n")
    targets = Tuple{String,String}[]
    for line in eachline(joinpath(ROOT, "Makefile"))
        m = match(r"^([a-zA-Z_-]+):.*?## (.*)$", line)
        m === nothing || push!(targets, (m[1], m[2]))
    end
    for (target, description) in sort(targets)
        printstyled(rpad(target, 20); color=:cyan)
        println(" ", description)
    end
    return nothing
end

"Install the development tools and instantiate the package and examples environments."
function install()
    Pkg.activate(DEV_ENV; shared=true)
    Pkg.add(DEV_TOOLS)
    for project in (ROOT, joinpath(ROOT, "examples"))
        Pkg.activate(project)
        Pkg.instantiate()
        Pkg.precompile()
    end
    return nothing
end

"""
Run all examples (examples/*.jl), each in its own Julia process with the examples
environment. The environment is resolved first, so that changes of the dependencies of
EKRem are taken into account.
"""
function examples()
    dir = joinpath(ROOT, "examples")
    Pkg.activate(dir)
    Pkg.resolve()
    Pkg.instantiate()
    scripts = sort(filter(f -> endswith(f, ".jl"), readdir(dir)))
    failed = String[]
    for script in scripts
        printstyled("\n▶ ", script, "\n"; color=:cyan, bold=true)
        cmd = `$(Base.julia_cmd()) --project=$dir $(joinpath(dir, script))`
        success(pipeline(cmd; stdout=stdout, stderr=stderr)) || push!(failed, script)
    end
    isempty(failed) || (println(stderr, "\nFailed examples: ", join(failed, ", ")); exit(1))
    println("\nAll $(length(scripts)) examples completed; figures are in examples/output/.")
    return nothing
end

"Run the tests with coverage tracking, print the coverage of src/ and write lcov.info."
function coverage()
    Pkg.test("EKRem"; coverage=true)
    files = Coverage.process_folder(joinpath(ROOT, "src"))
    println("\nLine coverage of src/:")
    for fc in sort(files; by=f -> f.filename)
        covered, total = Coverage.get_summary(fc)
        total == 0 && continue
        name = relpath(fc.filename, ROOT)
        println(
            "  ",
            rpad(name, 40),
            lpad(covered, 5),
            " / ",
            lpad(total, 5),
            "  (",
            round(100 * covered / total; digits=1),
            " %)",
        )
    end
    covered, total = Coverage.get_summary(files)
    println(
        "  ",
        rpad("Total", 40),
        lpad(covered, 5),
        " / ",
        lpad(total, 5),
        "  (",
        round(100 * covered / total; digits=1),
        " %)",
    )
    Coverage.LCOV.writefile(joinpath(ROOT, "lcov.info"), files)
    println("\nLCOV report written to lcov.info. Remove the coverage files with `make clean`.")
    return nothing
end

"Check the package quality with Aqua.jl (ambiguities, exports, dependencies, compat, piracy)."
function lint()
    Aqua.test_all(EKRem)
    return nothing
end

"Format the code, or only check the formatting if `check` is true."
function format_code(; check::Bool)
    formatted = JuliaFormatter.format(code_dirs(); overwrite=!check, verbose=check)
    if check && !formatted
        println(stderr, "Some files are not formatted. Run `make format`.")
        exit(1)
    end
    println(formatted ? "All files are formatted." : "Files have been formatted.")
    return nothing
end

"Increase the patch, minor or major version number in Project.toml."
function bump(level::AbstractString)
    file = joinpath(ROOT, "Project.toml")
    text = read(file, String)
    m = match(r"(?m)^version = \"(.+)\"$", text)
    m === nothing && error("No version entry found in Project.toml")
    v = VersionNumber(m[1])
    new = if level == "patch"
        VersionNumber(v.major, v.minor, v.patch + 1)
    elseif level == "minor"
        VersionNumber(v.major, v.minor + 1, 0)
    elseif level == "major"
        VersionNumber(v.major + 1, 0, 0)
    else
        error("Unknown version level \"$level\"; use patch, minor or major")
    end
    write(file, replace(text, m.match => "version = \"$new\""; count=1))
    println("Version bumped: $v → $new")
    return nothing
end

"Remove coverage (*.cov, lcov.info) and allocation (*.mem) files."
function clean()
    removed = 0
    for dir in code_dirs(), (root, _, files) in walkdir(dir), f in files
        if endswith(f, ".cov") || endswith(f, ".mem")
            rm(joinpath(root, f))
            removed += 1
        end
    end
    lcov = joinpath(ROOT, "lcov.info")
    isfile(lcov) && (rm(lcov); removed += 1)
    println("Removed $removed file(s).")
    return nothing
end

if TASK == "help"
    help()
elseif TASK == "install"
    install()
elseif TASK == "examples"
    examples()
elseif TASK == "coverage"
    coverage()
elseif TASK == "lint"
    lint()
elseif TASK == "format"
    format_code(; check=false)
elseif TASK == "format-check"
    format_code(; check=true)
elseif TASK == "bump"
    bump(TASK_ARG)
elseif TASK == "clean"
    clean()
else
    println(stderr, "Unknown task \"$TASK\". Run `make help`.")
    exit(1)
end
