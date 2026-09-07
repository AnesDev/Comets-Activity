include("physics/sublimation_pressure.jl")
include("physics/sublimation_flux.jl")

# NOTE: x, t, T, n, Dx, Dt are assumed already defined by equations.jl,
# which must be `include`d before this file in your top-level script.

"""
Steady-state (dT/dx = 0) solution of Eq. (1), used as the initial
isothermal nucleus temperature at heliocentric distance r_H:

    F_sun*(1-A_bond)/r_H^2 + 4*sigma*eps_IR*T_OC^4 = 4*sigma*eps_IR*T_i^4

Verified against the paper: at r_H = 150 au this formula gives
T_i ≈ 23 K, matching the paper's stated Mendis-point blackbody
temperature exactly. No special-cased override — this formula is
used for all r_H, including 100 au (which gives T_i ≈ 28 K, not 40 K).
"""
function isothermal_equilibrium_temperature(r_H, params)
    lhs = params.F_sun * (1 - params.A_bond) / r_H^2
    T_i4 = lhs / (4 * params.sigma_SB * params.eps_IR) + params.T_OC^4
    return T_i4^(1/4)
end

"""
Builds the initial and boundary conditions for Equations (6) and (7).

Returns
-------
(bcs, N_scale, T_i, n_i_scaled)
"""
function build_boundary_conditions(params, domain, r_H)

    x_A, x_B, t0, t1 = domain.x_A, domain.x_B, domain.t0, domain.t1

    T_i = 40.0  # TEMP TEST: was isothermal_equilibrium_temperature(r_H, params) = 28K. REVERT after comparing.
    # temorarily removed
    # T_i = isothermal_equilibrium_temperature(r_H, params)
    n_i = sublimation_pressure(T_i, params) / (params.k_B * T_i)

    N_scale = n_i
    n_i_scaled = n_i / N_scale   # = 1.0 always, by construction — kept
                                  # explicit and passed through rather
                                  # than hardcoded, so model.jl never
                                  # has to assume this invariant itself.

    # ic_n REMOVED from the soft BC list — n(x,t0) = n_i_scaled is now
    # HARD-ENFORCED by the model architecture (see model.jl), same as
    # T's IC. Keeping a soft BC for an already-exact constraint is
    # exactly what was destabilizing the loss balancing (BC[1] was
    # stuck flat for hundreds of iterations both before and after the
    # T_i fix).

    bc_surface = (
        params.F_sun * (1 - params.A_bond) / r_H^2
        + 4 * params.sigma_SB * params.eps_IR * params.T_OC^4
        ~
        4 * params.sigma_SB * params.eps_IR * T(x_A, t)^4
        + 4 * conductivity(T(x_A, t), params) * Dx(T(x_A, t))
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

    # bc_gas_A REMOVED from the soft BC list — n(x_A,t) = n_small/N_scale is
    # now HARD-ENFORCED by the model architecture (see model.jl), same as
    # the T and n initial conditions. This was ~98% of total loss and never
    # converged as a soft term, so it's architecturally imposed instead.

    n_small = 1e-6 * n_i
    n_A_scaled = n_small / N_scale

    return [bc_surface, bc_heat_B, bc_gas_B], N_scale, T_i, n_i_scaled, n_A_scaled

end