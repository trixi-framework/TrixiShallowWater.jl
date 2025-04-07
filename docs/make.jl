using TrixiShallowWater
using Documenter
using DocumenterInterLinks
using Literate

# Provide external links to the Trixi.jl docs (project root and inventory file)
links = InterLinks("Trixi" => ("https://trixi-framework.github.io/Trixi.jl/stable/",
                               "https://trixi-framework.github.io/Trixi.jl/stable/objects.inv"))

# Create tutorial section with Literature
TUTORIAL_DIR = joinpath(@__DIR__, "src", "tutorials")
OUTPUT_DIR = joinpath(@__DIR__, "src", "tutorials")

tutorial_list = [
    "elixir_shallowwater_dam_break_triangular.jl"
]

tutorial_pages = [
    "Dam Break" => "tutorials/elixir_shallowwater_dam_break_triangular.md"
]

# Create markdown files
for tutorial in tutorial_list
    Literate.markdown(joinpath(TUTORIAL_DIR, tutorial), OUTPUT_DIR;)
end

# Copy list of authors to not need to synchronize it manually
authors_text = read(joinpath(dirname(@__DIR__), "AUTHORS.md"), String)
authors_text = replace(authors_text,
                       "in the [LICENSE.md](LICENSE.md) file" => "under [License](@ref)")
write(joinpath(@__DIR__, "src", "authors.md"), authors_text)

# Copy contributing information to not need to synchronize it manually
contributing_text = read(joinpath(dirname(@__DIR__), "CONTRIBUTING.md"), String)
contributing_text = replace(contributing_text,
                            "[LICENSE.md](LICENSE.md)" => "[License](@ref)",
                            "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref)")
write(joinpath(@__DIR__, "src", "contributing.md"), contributing_text)

# Copy code of conduct to not need to synchronize it manually
code_of_conduct_text = read(joinpath(dirname(@__DIR__), "CODE_OF_CONDUCT.md"), String)
code_of_conduct_text = replace(code_of_conduct_text,
                                "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref)")
write(joinpath(@__DIR__, "src", "code_of_conduct.md"), code_of_conduct_text)

# Copy contents form README to the starting page to not need to synchronize it manually
readme_text = read(joinpath(dirname(@__DIR__), "README.md"), String)
readme_text = replace(readme_text,
                       "[LICENSE.md](LICENSE.md)" => "[License](@ref)",
                       "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref)",
                       "<p" => "```@raw html\n<p",
                       "p>" => "p>\n```",
                       r"\[comment\].*\n" => "")    # remove comments
write(joinpath(@__DIR__, "src", "home.md"), readme_text)

DocMeta.setdocmeta!(TrixiShallowWater, :DocTestSetup, :(using TrixiShallowWater);
                    recursive = true)

# TODO: create changelog

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
         pages = ["Home" => "home.md",
             "Installation" => "installation.md",
             "Tutorials" => tutorial_pages,
             "Advanced topics & developers" => ["Development" => "development.md"
                                                "Testing" => "testing.md"],
             "Authors" => "authors.md",
             "Contributing" => "contributing.md",
             "Code of Conduct" => "code_of_conduct.md",
             "License" => "license.md",
             "Index" => "index.md"],
         plugins = [links],)

deploydocs(repo = "github.com/trixi-framework/TrixiShallowWater.jl",
           devbranch = "main",
           push_preview = true)
