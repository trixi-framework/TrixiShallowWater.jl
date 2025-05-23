using TrixiShallowWater
using Documenter
using DocumenterInterLinks
using Literate
using Changelog: Changelog

# Provide external links to the Trixi.jl docs (project root and inventory file)
links = InterLinks("Trixi" => ("https://trixi-framework.github.io/TrixiDocumentation/stable/",
                               "https://trixi-framework.github.io/TrixiDocumentation/stable/objects.inv"))

# Create tutorial section with Literature
TUTORIAL_DIR = joinpath(@__DIR__, "src", "tutorials")
OUTPUT_DIR = joinpath(@__DIR__, "src", "tutorials")

tutorial_list = [
    "elixir_shallowwater_dam_break_triangular.jl",
    "elixir_shallowwater_monai_tsunami.jl"
]

tutorial_pages = [
    "Introduction" => "tutorials/introduction.md",
    "Dam break over triangular bottom topograhy" => "tutorials/elixir_shallowwater_dam_break_triangular.md",
    "Okushiri tsunami" => "tutorials/elixir_shallowwater_monai_tsunami.md"
]

# Create markdown files
for tutorial in tutorial_list
    Literate.markdown(joinpath(TUTORIAL_DIR, tutorial), OUTPUT_DIR;)
end

# Copy list of authors to not need to synchronize it manually.
# Since the authors header exists twice we create a unique identifier for the docs section.
authors_text = read(joinpath(dirname(@__DIR__), "AUTHORS.md"), String)
authors_text = replace(authors_text,
                       "in the [LICENSE.md](LICENSE.md) file" => "under [License](@ref)",
                       "# Authors" => "# [Authors](@id trixi_sw_authors)")
write(joinpath(@__DIR__, "src", "authors.md"), authors_text)

# Copy contributing information to not need to synchronize it manually
contributing_text = read(joinpath(dirname(@__DIR__), "CONTRIBUTING.md"), String)
contributing_text = replace(contributing_text,
                            "[LICENSE.md](LICENSE.md)" => "[License](@ref)",
                            "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref trixi_sw_authors)")
write(joinpath(@__DIR__, "src", "contributing.md"), contributing_text)

# Copy code of conduct to not need to synchronize it manually
code_of_conduct_text = read(joinpath(dirname(@__DIR__), "CODE_OF_CONDUCT.md"), String)
code_of_conduct_text = replace(code_of_conduct_text,
                               "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref trixi_sw_authors)")
write(joinpath(@__DIR__, "src", "code_of_conduct.md"), code_of_conduct_text)

# Copy contents form README to the starting page to not need to synchronize it manually
readme_text = read(joinpath(dirname(@__DIR__), "README.md"), String)
readme_text = replace(readme_text,
                      "[LICENSE.md](LICENSE.md)" => "[License](@ref)",
                      "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref trixi_sw_authors)",
                      "<p" => "```@raw html\n<p",
                      "p>" => "p>\n```",
                      r"\[comment\].*\n" => "")    # remove comments
write(joinpath(@__DIR__, "src", "home.md"), readme_text)

DocMeta.setdocmeta!(TrixiShallowWater, :DocTestSetup, :(using TrixiShallowWater);
                    recursive = true)

# Create changelog
Changelog.generate(Changelog.Documenter(),                        # output type
                   joinpath(@__DIR__, "..", "NEWS.md"),           # input file
                   joinpath(@__DIR__, "src", "changelog_tmp.md"); # output file
                   repo = "trixi-framework/TrixiShallowWater.jl", # default repository for links
                   branch = "main",)
# Fix edit URL of changelog
open(joinpath(@__DIR__, "src", "changelog.md"), "w") do io
    for line in eachline(joinpath(@__DIR__, "src", "changelog_tmp.md"))
        if startswith(line, "EditURL")
            line = "EditURL = \"https://github.com/trixi-framework/TrixiShallowWater.jl/blob/main/NEWS.md\""
        end
        println(io, line)
    end
end
# Remove temporary file
rm(joinpath(@__DIR__, "src", "changelog_tmp.md"))

# Make documentation
makedocs(;
         modules = [TrixiShallowWater],
         authors = "Andrew R. Winters <andrew.ross.winters@liu.se>, Michael Schlottke-Lakemper <michael@sloede.com>",
         repo = "https://github.com/trixi-framework/TrixiShallowWater.jl/blob/{commit}{path}#{line}",
         sitename = "TrixiShallowWater.jl",
         format = Documenter.HTML(;
                                  # Disable pretty URLs during manual testing
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  # Set canonical URL to GitHub pages URL
                                  canonical = "https://trixi-framework.github.io/TrixiShallowWater.jl",
                                  edit_link = "main",
                                  size_threshold_ignore = ["reference.md"],),
         # Explicitly specify documentation structure
         pages = ["Home" => "index.md",
             "Installation" => "installation.md",
             "Tutorials" => tutorial_pages,
             "Advanced topics & developers" => ["Development" => "development.md",
                 "Style guide" => "styleguide.md",
                 "Testing" => "testing.md"],
             "Changelog" => "changelog.md",
             "Authors" => "authors.md",
             "Contributing" => "contributing.md",
             "Code of Conduct" => "code_of_conduct.md",
             "License" => "license.md",
             "Reference" => "reference.md"],
         plugins = [links],)

deploydocs(repo = "github.com/trixi-framework/TrixiShallowWater.jl",
           devbranch = "main",
           push_preview = true)
