using TrixiShallowWater
using Documenter
using DocumenterInterLinks

# Provide external links to the Trixi.jl docs (project root and inventory file)
links = InterLinks("Trixi" => ("https://trixi-framework.github.io/Trixi.jl/stable/",
                               "https://trixi-framework.github.io/Trixi.jl/stable/objects.inv"))

# Copy list of authors to not need to synchronize it manually
authors_text = read(joinpath(dirname(@__DIR__), "AUTHORS.md"), String)
authors_text = replace(authors_text,
                       "in the [LICENSE.md](LICENSE.md) file" => "under [License](@ref)")
write(joinpath(@__DIR__, "src", "authors.md"), authors_text)

DocMeta.setdocmeta!(TrixiShallowWater, :DocTestSetup, :(using TrixiShallowWater);
                    recursive = true)

# Copy some files from the repository root directory to the docs and modify them
# as necessary
# Based on: https://github.com/ranocha/SummationByPartsOperators.jl/blob/0206a74140d5c6eb9921ca5021cb7bf2da1a306d/docs/make.jl#L27-L41
open(joinpath(@__DIR__, "src", "code_of_conduct.md"), "w") do io
    # Point to source license file
    println(io,
            """
            ```@meta
            EditURL = "https://github.com/trixi-framework/TrixiShallowWater.jl/blob/main/CODE_OF_CONDUCT.md"
            ```
            """)
    # Write the modified contents
    println(io, "# [Code of Conduct](@id code-of-conduct)")
    println(io, "")
    for line in eachline(joinpath(dirname(@__DIR__), "CODE_OF_CONDUCT.md"))
        line = replace(line, "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref)")
        println(io, "> ", line)
    end
end

open(joinpath(@__DIR__, "src", "contributing.md"), "w") do io
    # Point to source license file
    println(io,
            """
            ```@meta
            EditURL = "https://github.com/trixi-framework/TrixiTrixiShallowWater.jl/blob/main/CONTRIBUTING.md"
            ```
            """)
    # Write the modified contents
    for line in eachline(joinpath(dirname(@__DIR__), "CONTRIBUTING.md"))
        line = replace(line, "[LICENSE.md](LICENSE.md)" => "[License](@ref)")
        line = replace(line, "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref)")
        println(io, line)
    end
end

makedocs(;
         modules = [TrixiShallowWater],
         authors = "Andrew R. Winters <andrew.ross.winters@liu.se>, Michael Schlottke-Lakemper <michael@sloede.com>",
         repo = "https://github.com/trixi-framework/TrixiShallowWater.jl/blob/{commit}{path}#{line}",
         sitename = "TrixiShallowWater.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://trixi-framework.github.io/TrixiShallowWater.jl",
                                  edit_link = "main",
                                  assets = String[],),
         pages = [
             "Home" => "index.md",
             "Installation" => "installation.md",
             "Tutorials" => "tutorial.md",
             "Advanced topics & developers" => ["Development" => "development.md"
                                                "Testing" => "testing.md"],
             "Authors" => "authors.md",
             "Contributing" => "contributing.md",
             "Code of Conduct" => "code_of_conduct.md",
             "License" => "license.md",
         ],
         plugins = [links],)

deploydocs(;
           repo = "github.com/trixi-framework/TrixiShallowWater.jl",
           devbranch = "main",)
