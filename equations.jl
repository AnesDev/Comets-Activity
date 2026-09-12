using ModelingToolkit

include("physics/pressure.jl")
include("physics/heat_capacity.jl")
include("physics/conductivity.jl")
include("physics/thermal_speed.jl")
include("physics/viscosity.jl")
include("physics/permeability.jl")
include("physics/knudsen_diffusion.jl")
include("physics/sublimation_flux.jl")

@parameters x t

@variables T(..) n(..)

Dt = Differential(t)
Dx = Differential(x)

function build_equations(params, N_scale, T_i)

    heat_eq_scale = params.latent_heat_CO * sublimation_flux(T_i, params)

    eq_heat = (
        params.density *
        heat_capacity(T(x, t)) *
        Dt(T(x, t)) / heat_eq_scale
        ~
        Dx(
            conductivity(T(x, t), params) *
            Dx(T(x, t))
        ) / heat_eq_scale
    )

    dn_dx = Dx(n(x, t))
    dT_dx = Dx(T(x, t))

    n_phys     = N_scale * n(x, t)
    dn_phys_dx = N_scale * dn_dx

    gas_flux_scale = sublimation_flux(T_i, params) / params.CO_molecular_mass
    gas_rate_scale = gas_flux_scale / params.mantle_thickness

    eq_gas = (
        params.porosity * Dt(n(x, t)) * N_scale / gas_rate_scale
        ~
        Dx(
            knudsen_diffusion(T(x, t), params) * dn_phys_dx
            +
            permeability(params) / viscosity(T(x, t), params) *
            n_phys *
            (
                params.k_B * (dn_phys_dx * T(x, t) + n_phys * dT_dx)
            )
        ) / gas_rate_scale
    )

    return eq_heat, eq_gas, heat_eq_scale, gas_rate_scale

end