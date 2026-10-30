# =============================================================================
#  Grid
# =============================================================================

const _TICK_BASE = IdDict{Any,String}()

"""
    _tick!(bar, vals)

Advance a progress bar by one step, keeping `vals` visible.

`ProgressMeter`'s `showvalues` puts the values on extra lines *below* the bar and
redraws them with cursor-up escapes.  Jupyter honours a carriage return but not
cursor movement, so under IJulia those lines pile up instead of being overwritten
-- a fresh block on every step rather than one updating bar.  A single-line bar
does render correctly there (`thermodynamics!` uses one and is fine, precisely
because it passes no `showvalues`), so in a notebook the values are folded into
the bar's own description instead of being dropped.  In a terminal `showvalues`
is used as intended.

Note the bars here differ in *granularity*, which is not a rendering problem:
`qnms` ticks once per grid in `Ns`, so three ticks and it is done, while
`thermodynamics!` ticks `nphi` times and looks smoothly animated.
"""
function _tick!(bar, vals)
    isdefined(Main, :IJulia) || return next!(bar; showvalues = vals)
    base = get!(_TICK_BASE, bar, bar.desc)
    return next!(bar; desc = string(base, "(",
                                    join([string(k, " ", v) for (k, v) in vals], ", "),
                                    ") "))
end

"""
    cheby(N) -> x, D, D2

`N` Chebyshev-Gauss-Lobatto points mapped to [0,1] (`x[1] = 0` boundary,
`x[N] = 1` horizon) with the first and second derivative matrices.
"""
function cheby(N::Int)
    x = -cos.(π * (0:N-1) / (N - 1))
    c = [i == 0 || i == N-1 ? 2.0 : 1.0 for i = 0:N-1]
    c = c .* (-1).^(0:N-1)
    X = repeat(x, 1, N); dX = X - X'
    D = (c ./ c') ./ (dX + I); D = D - Diagonal(vec(sum(D, dims = 2)))
    x = (x .+ 1) / 2; D = 2 * D; D2 = D * D
    return x, D, D2
end
