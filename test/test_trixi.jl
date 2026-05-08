using TrixiTest

macro test_trixi_include(expr, args...)
    local add_to_additional_ignore_content = [
        # We need to ignore steady state information reported by our callbacks
        r"┌ Info:   Steady state tolerance reached\n│   steady_state_callback .+\n└   t = .+\n",
        # NOTE: These warnings arose from Julia 1.10 onwards
        r"WARNING: Method definition .* in module .* at .* overwritten .*.\n"
    ]
    args = append_to_kwargs(args, :additional_ignore_content,
                            add_to_additional_ignore_content)
    ex = quote
        @test_trixi_include_base($expr, $(args...))
    end
    return esc(ex)
end
