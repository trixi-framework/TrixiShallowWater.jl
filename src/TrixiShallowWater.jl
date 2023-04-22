module TrixiShallowWater

# import @reexport now to make it available for further imports/exports
using Reexport: @reexport

# Make all of Trixi.jl available to a user of this package
@reexport using Trixi

# Write your package code here.
foo() = true
bar() = false

end
