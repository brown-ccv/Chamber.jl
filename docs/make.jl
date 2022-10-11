push!(LOAD_PATH, "../src/")

using Documenter
using Chamber

makedocs(;
    sitename="Chamber.jl",
    format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules=[Chamber],
    pages = [
        "Introduction" => "index.md",
        "Background" => "background.md",
        # Add more pages as needed here
    ],
)

# Hosting: update repo
deploydocs(;
    repo="github.com/CallieHsu/Chamber.jl.git",
    push_preview=true,
)
