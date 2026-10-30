# =============================================================================
#  Containers
# =============================================================================

"""One background solution on one grid, plus its thermodynamic data."""
struct Background
    N::Int
    x::Vector{Float64}
    D::Matrix{Float64}
    D2::Matrix{Float64}
    f0::Vector{Float64}
    g0::Vector{Float64}
    phi0::Vector{Float64}
    phiH::Float64
    T::Float64
    S::Float64
    Lambda::Float64
    res::Float64
    converged::Bool
end

"""
Outcome of a QNM computation at one (phiH, q).

`raw[i]` is the full unfiltered spectrum on grid `Ns[i]`; `omega`/`lambda` are
the frequencies that survived the multi-resolution test, `err` their drift
between the two finest grids (so `-log10(err/|lambda|)` ~ significant digits).
"""
struct QNMResult
    phiH::Float64
    q::Float64
    T::Float64
    Ns::Vector{Int}
    raw::Vector{Vector{ComplexF64}}
    omega::Vector{ComplexF64}
    lambda::Vector{ComplexF64}
    err::Vector{Float64}
end

"""
    Bea_Mateos(phiM, phiQ; chi = 1.0, nu = 1.0)

The Bea-Mateos model with potential parameters `phiM`, `phiQ`.  The constructor
builds V, V', V'' and locates `phiMax`, the minimum of V, which bounds the
physical range of the horizon value phiH.

Fields worth reading: `V`, `Vp`, `Vpp`, `Vfuncs`, `phiMax`, `thermo`, `Tc`.
Solved backgrounds are cached in `bg`, keyed by `(N, phiH)`.
"""
mutable struct Bea_Mateos
    phiM::Float64
    phiQ::Float64
    chi::Float64
    nu::Float64

    V::Function
    Vp::Function
    Vpp::Function
    Vfuncs::Vector{Function}
    phiMax::Float64

    bg::Dict{Tuple{Int,Float64},Background}
    thermo::Union{Nothing,NamedTuple}
    Tc::Union{Nothing,Float64}
    qnm_runs::Vector{QNMResult}
end

function Bea_Mateos(phiM::Real, phiQ::Real; chi::Real = 1.0, nu::Real = 1.0,
                    phisearch::Real = 5.0)
    M, Q = float(phiM), float(phiQ)

    V(p)   = -1.5p^2 - (1/3)p^4 + p^6*(1/(2M^4) - 1/(3M^2) - 2/Q) +
             p^8*(-(72M^2 - 16M^4 + Q)/(12M^4*Q)) +
             p^10*((2*(27 + Q/M^2))/(3Q^2)) - (4/3)p^12/Q^2
    Vp(p)  = -3p - (4/3)p^3 + 6p^5*(1/(2M^4) - 1/(3M^2) - 2/Q) +
             8p^7*(-(72M^2 - 16M^4 + Q)/(12M^4*Q)) +
             10p^9*(2*(27 + Q/M^2)/(3Q^2)) - 16p^11/Q^2
    Vpp(p) = -3 - 4p^2 + 30p^4*(1/(2M^4) - 1/(3M^2) - 2/Q) +
             56p^6*(-(72M^2 - 16M^4 + Q)/(12M^4*Q)) +
             90p^8*(2*(27 + Q/M^2)/(3Q^2)) - 176p^10/Q^2

    phiMax = Optim.minimizer(optimize(V, 0.0, float(phisearch)))
    if phiMax > 0.99 * phisearch
        @warn "phiMax hit the edge of the search window; raise `phisearch`." phiMax phisearch
    end

    return Bea_Mateos(M, Q, float(chi), float(nu), V, Vp, Vpp,
                      Function[V, Vp, Vpp], phiMax,
                      Dict{Tuple{Int,Float64},Background}(), nothing, nothing,
                      QNMResult[])
end

function Base.show(io::IO, m::Bea_Mateos)
    print(io, "Bea_Mateos(phiM = ", m.phiM, ", phiQ = ", m.phiQ,
          ", chi = ", m.chi, ", nu = ", m.nu, ")\n")
    print(io, "  phiMax = ", round(m.phiMax, digits = 6),
          "   V(phiMax) = ", round(m.V(m.phiMax), digits = 6), "\n")
    print(io, "  cached backgrounds: ", length(m.bg),
          "   thermodynamics: ", m.thermo === nothing ? "not computed" :
              string(length(m.thermo.phiH), " points"),
          "   Tc: ", m.Tc === nothing ? "-" : string(round(m.Tc, digits = 6)), "\n")
    print(io, "  QNM runs: ", length(m.qnm_runs))
end
