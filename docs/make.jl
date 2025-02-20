using TrixiShallowWater
using Documenter
using DocumenterInterLinks

# Provide external links to the Trixi.jl docs (project root and inventory file)
links = InterLinks("Trixi" => ("https://trixi-framework.github.io/Trixi.jl/stable/",
                               "https://trixi-framework.github.io/Trixi.jl/stable/objects.inv"))

DocMeta.setdocmeta!(TrixiShallowWater, :DocTestSetup, :(using TrixiShallowWater);
                    recursive = true)

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
         pages = ["Home" => "index.md"],
         plugins = [links],)

deploydocs(;
           repo = "github.com/trixi-framework/TrixiShallowWater.jl",
           devbranch = "main",)
