using ModelingToolkit

include("physics/pressure.jl")
include("physics/heat_capacity.jl")
include("physics/conductivity.jl")
include("physics/thermal_speed.jl")
include("physics/viscosity.jl")
include("physics/permeability.jl")
include("physics/knudsen_diffusion.jl")

@parameters x t

@variables T(..) n(..)   # n(x,t) is now DIMENSIONLESS: n(x,t) ≈ O(1).
                          # Physical density = N_scale * n(x,t).

Dt = Differential(t)
Dx = Differential(x)

function build_equations(params, N_scale)

    eq_heat = (
        params.density *
        heat_capacity(T(x, t)) *
        Dt(T(x, t))
        ~
        Dx(
            conductivity(T(x, t), params) *
            Dx(T(x, t))
        )
    )

    dn_dx = Dx(n(x, t))            # Dx applied to BARE n(x,t) only
    dT_dx = Dx(T(x, t))            # Dx applied to BARE T(x,t) only

    n_phys     = N_scale * n(x, t)
    dn_phys_dx = N_scale * dn_dx   # scale AFTER differentiating — never before

    eq_gas = (
        params.porosity * Dt(n(x, t))
        ~
        (1 / N_scale) * Dx(
            knudsen_diffusion(T(x, t), params) * dn_phys_dx
            +
            permeability(params) / viscosity(T(x, t), params) *
            n_phys *
            (
                params.k_B * (dn_phys_dx * T(x, t) + n_phys * dT_dx)
            )
        )
    )

    return eq_heat, eq_gas

end