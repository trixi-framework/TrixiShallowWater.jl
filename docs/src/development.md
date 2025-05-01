# Development

This page contains some helpful information for the development of TrixiShallowWater.jl. Further
information about helpful tools for package development in Julia can be found on the
[development page](https://trixi-framework.github.io/TrixiDocumentation/stable/development/) of the Trixi.jl docs.

## Releasing a new version of TrixiShallowWater

- Check whether everything is okay, tests pass etc.
- Set the new version number in `Project.toml` according to the Julian version of semver.
  Commit and push.
- Comment `@JuliaRegistrator register` on the commit setting the version number.
- `JuliaRegistrator` will create a PR with the new version in the General registry.
  Wait for it to be merged.
- Increment the version number in `Project.toml` again with suffix `-DEV`. For example,
  if you have released version `v0.2.0`, use `v0.2.1-DEV` as new version number.



## Preview the documentation

You can build the documentation of TrixiShallowWater.jl locally by running
```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs --color=yes docs/make.jl
```
from the TrixiShallowWater.jl main directory. Then, you can look at the html files generated in
`docs/build`.
For PRs triggered from branches inside the TrixiShallowWater.jl main repository previews of
the new documentation are generated at
`https://trixi-framework.github.io/TrixiShallowWater.jl/previews/PRXXX`,
where `XXX` is the number of the PR.
Note, this does not work for PRs from forks for security reasons (since anyone could otherwise push
arbitrary stuff, including malicious code).


## Developing with a local Trixi.jl version

TrixiShallowWater.jl has Trixi.jl as a dependency and uses Trixi.jl's implementation.
When developing TrixiShallowWater.jl, one may want to change functions in Trixi.jl to allow them to be reused
in TrixiShallowWater.jl.
To use a locally modified Trixi.jl clone instead of a Trixi.jl release, one can tell Pkg
to use the local source code of Trixi.jl instead of a registered version by running
```julia-repl
using Pkg
Pkg.develop(PackageSpec(path="path/to/Trixi.jl"))
```

To switch back from a local version to a Trixi.jl release run
```julia-repl
Pkg.free("Trixi")
```
