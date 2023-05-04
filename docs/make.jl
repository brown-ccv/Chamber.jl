push!(LOAD_PATH, "../src/")

using Documenter
using Chamber

makedocs(;
    sitename="Chamber.jl",
    format=Documenter.HTML(prettyurls=get(ENV, "CI", nothing) == "true"),
    modules=[Chamber]
)

# Hosting: update repo
deploydocs(;
    repo="github.com/brown-ccv/Chamber.jl.git",
    push_preview=true
)
