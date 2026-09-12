include("physics/sublimation_pressure.jl")
include("physics/sublimation_flux.jl")

# NOTE: x, t, T, n, Dx, Dt are assumed already defined by equations.jl,
# which must be `include`d before this file in your top-level script.

function isothermal_equilibrium_temperature(r_H, params)
    lhs = params.F_sun * (1 - params.A_bond) / r_H^2
    T_i4 = lhs / (4 * params.sigma_SB * params.eps_IR) + params.T_OC^4
    return T_i4^(1/4)
end

"""
Builds the initial and boundary conditions for Equations (6) and (7).
"""
function build_boundary_conditions(params, domain, r_H; T_i_override=nothing)

    x_A, x_B, t0, t1 = domain.x_A, domain.x_B, domain.t0, domain.t1

    T_i = T_i_override === nothing ? isothermal_equilibrium_temperature(r_H, params) : T_i_override

    n_i = sublimation_pressure(T_i, params) / (params.k_B * T_i)

    N_scale = n_i
    n_i_scaled = n_i / N_scale 

    solar_flux_scale = params.F_sun * (1 - params.A_bond) / r_H^2

    bc_surface = (
        (params.F_sun * (1 - params.A_bond) / r_H^2
        + 4 * params.sigma_SB * params.eps_IR * params.T_OC^4) / solar_flux_scale
        ~
        (4 * params.sigma_SB * params.eps_IR * T(x_A, t)^4
        + 4 * conductivity(T(x_A, t), params) * Dx(T(x_A, t))) / solar_flux_scale
    )

    heat_flux_scale = params.latent_heat_CO * sublimation_flux(T_i, params)
    bc_heat_B = (
        conductivity(T(x_B,t),params)*Dx(T(x_B,t)) / heat_flux_scale
        ~
        params.latent_heat_CO * sublimation_flux(T(x_B,t), params) / heat_flux_scale
    )

    dn_dx_B = Dx(n(x_B, t))
    dT_dx_B = Dx(T(x_B, t))

    gas_flux_scale = sublimation_flux(T_i, params) / params.CO_molecular_mass

    bc_gas_B = (
        (
            knudsen_diffusion(T(x_B, t), params) * N_scale * dn_dx_B
            +
            permeability(params) / viscosity(T(x_B, t), params)
            * N_scale * n(x_B, t) *
            (
                params.k_B * (N_scale * dn_dx_B * T(x_B, t) + N_scale * n(x_B, t) * dT_dx_B)
            )
        ) / gas_flux_scale
        ~
        sublimation_flux(T(x_B, t), params) / params.CO_molecular_mass / gas_flux_scale
    )


    n_small = 1e-6 * n_i
    n_A_scaled = n_small / N_scale

    return [bc_surface, bc_heat_B, bc_gas_B], N_scale, T_i, n_i_scaled, n_A_scaled, solar_flux_scale, heat_flux_scale, gas_flux_scale

end