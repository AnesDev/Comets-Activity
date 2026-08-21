include("physics/sublimation_pressure.jl")
include("physics/sublimation_flux.jl")

# NOTE: x, t, T, n, Dx, Dt are assumed already defined by equations.jl,
# which must be `include`d before this file in your top-level script.

"""
Steady-state (dT/dx = 0) solution of Eq. (1), used as the initial
isothermal nucleus temperature at heliocentric distance r_H:

    F_sun*(1-A_bond)/r_H^2 + 4*sigma*eps_IR*T_OC^4 = 4*sigma*eps_IR*T_i^4


"""
# Supervisor-confirmed values, not derived — see email thread.
# r_H = 100 au -> T_i ≈ 40 K (confirmed explicitly)
# r_H = 10 au, 4 au -> NOT YET CONFIRMED, still using derived formula below

function isothermal_equilibrium_temperature(r_H, params)
    if r_H ≈ 100.0
        return 40.0  # supervisor-confirmed value, overrides derivation below
    end
    lhs = params.F_sun * (1 - params.A_bond) / r_H^2
    T_i4 = lhs / (4 * params.sigma_SB * params.eps_IR) + params.T_OC^4
    return T_i4^(1/4)
end

"""
Builds the initial and boundary conditions for Equations (6) and (7).

Arguments
---------
params : Params structure (see below for fields this function needs
         that must be added to Params — they aren't there yet).
domain : NamedTuple with x_A, x_B (spatial bounds, meters) and
         t0, t1 (time bounds, seconds).
r_H    : heliocentric distance in au. Treated as CONSTANT for now
         (per your supervisor's 3-case plan: 100, 10, 4 au). Once
         the Kepler solver exists, replace this with r_H(t) and turn
         the surface BC below into a genuinely time-dependent Robin
         condition instead of one evaluated at fixed r_H.

Returns
-------
Vector of ModelingToolkit equations (ICs + BCs), ready to hand to
NeuralPDE.jl's PDESystem alongside build_equations(params).
"""


function build_boundary_conditions(params, domain, r_H)

    x_A, x_B, t0, t1 = domain.x_A, domain.x_B, domain.t0, domain.t1

    T_i = isothermal_equilibrium_temperature(r_H, params)
    n_i = sublimation_pressure(T_i, params) / (params.k_B * T_i)

    N_scale = n_i   # normalize so the network's n-channel stays O(1)

    ic_T = T(x, t0) ~ T_i
    ic_n = n(x, t0) ~ n_i / N_scale   # = 1, exactly

    bc_surface = (
        params.F_sun * (1 - params.A_bond) / r_H^2
        + 4 * params.sigma_SB * params.eps_IR * params.T_OC^4
        ~
        4 * params.sigma_SB * params.eps_IR * T(x_A, t)^4
        + 4 * conductivity(T(x_A, t), params) * Dx(T(x_A, t))
    )

    bc_heat_B = (
            conductivity(T(x_B, t), params) * Dx(T(x_B, t))
            ~
            params.latent_heat_CO * sublimation_flux(T(x_B, t), params)
        )

    dn_dx_B = Dx(n(x_B, t))            # Dx applied to BARE n(x,t) only
    dT_dx_B = Dx(T(x_B, t))            # Dx applied to BARE T(x,t) only

    n_phys_B     = N_scale * n(x_B, t)
    dn_phys_dx_B = N_scale * dn_dx_B   # scale AFTER differentiating

    gas_flux_B = (
        knudsen_diffusion(T(x_B, t), params) * dn_phys_dx_B
        +
        permeability(params) / viscosity(T(x_B, t), params)
        * n_phys_B *
        (
            params.k_B * (dn_phys_dx_B * T(x_B, t) + n_phys_B * dT_dx_B)
        )
    )
    bc_gas_B = (
        (1 / N_scale) * gas_flux_B
        ~
        (1 / N_scale) * sublimation_flux(T(x_B, t), params) / params.CO_molecular_mass
    )


    n_small = 1e-6 * n_i
    bc_gas_A = n(x_A, t) ~ n_small / N_scale

    return [ic_T, ic_n, bc_surface, bc_heat_B, bc_gas_B, bc_gas_A], N_scale

end

#=
function build_boundary_conditions(params, domain, r_H)

    x_A, x_B, t0, t1 = domain.x_A, domain.x_B, domain.t0, domain.t1

    # ---- Initial conditions ----------------------------------------
    T_i = isothermal_equilibrium_temperature(r_H, params)
    n_i = sublimation_pressure(T_i, params) / (params.k_B * T_i)  # ideal gas, inverted

    ic_T = T(x, t0) ~ T_i
    ic_n = n(x, t0) ~ n_i

    # ---- Surface BC at x = A: Robin condition, Eq. (1) --------------
    # F_sun*(1-A_bond)/r_H^2 + 4*sigma*eps_IR*T_OC^4
    #     = 4*sigma*eps_IR*T(A,t)^4 + 4*k(T(A,t))*dT/dx|_A
    bc_surface = (
        params.F_sun * (1 - params.A_bond) / r_H^2
        + 4 * params.sigma_SB * params.eps_IR * params.T_OC^4
        ~
        4 * params.sigma_SB * params.eps_IR * T(x_A, t)^4
        + 4 * conductivity(T(x_A, t), params) * Dx(T(x_A, t))
    )
    # NOTE: this equation should hold for all t in [t0, t1] at x = x_A.
    # Once r_H becomes r_H(t), swap the constant r_H^2 above for r_H(t)^2.

    # ---- Ice-front BC at x = B: Neumann conditions ------------------

    # Heat, Eq. (3): conducted flux == energy consumed by sublimation
    #
    # UNIT FIX: params.E_CO (7450 J/mol) is the Arrhenius activation
    # energy from Eq. (5) — it is NOT the same quantity as the latent
    # heat needed here. The Appendix confirms this explicitly: "L =
    # 2.66e5 J/kg is the latent heat of CO from Table 1" — that's
    # E_CO converted from J/mol to J/kg via the molar mass, not the
    # raw 7450 value. Eq. (3)'s dimensional balance only works with
    # J/kg: k(T)*dT/dx [W/m^2] = L [J/kg] * f_CO [kg/(m^2 s)].
    # See params.latent_heat_CO (added to Params below).
    bc_heat_B = (
        conductivity(T(x_B, t), params) * Dx(T(x_B, t))
        ~
        params.latent_heat_CO * sublimation_flux(T(x_B, t), params)
    )

    # Gas, Eq. (4)/(10): net number flux at the ice front is bounded
    # above by f_CO. f_CO is a MASS flux (kg/m^2/s); Eq. (7) is written
    # in terms of number density n, so convert via CO_molecular_mass.
    #
    # Sign convention: this uses the SAME (unsigned) flux expression
    # that appears inside Dx(...) in eq_gas in equations.jl. That's
    # the thing that actually matters for a self-consistent PINN — the
    # interior PDE residual and this boundary residual must reference
    # the same physical quantity. Do not flip the sign here alone; if
    # a sign correction is ever needed it has to be applied to BOTH
    # eq_gas and this BC together, after re-deriving the direction
    # convention from Eq. (10) — the paper's own text doesn't show an
    # explicit minus sign either.
    gas_flux_B = (
        knudsen_diffusion(T(x_B, t), params) * Dx(n(x_B, t))
        +
        permeability(params) / viscosity(T(x_B, t), params)
        * n(x_B, t) * 
        (
            params.k_B * (Dx(n(x_B, t)) * T(x_B, t) + n(x_B, t) * Dx(T(x_B, t)))
        )
    )

    bc_gas_B = (
        gas_flux_B ~ sublimation_flux(T(x_B, t), params) / params.CO_molecular_mass
    )

    # ---- Surface BC at x = A for the GAS equation: Dirichlet ---------
    # Confirmed with supervisor: no derived physical condition here —
    # the surface is (approximately) open to vacuum, and the model is
    # not sensitive to the exact value as long as it's small. So this
    # is Dirichlet, not Neumann, fixed to n_small.
    #
    # n_small is chosen as the equilibrium density at the (cold) surface
    # temperature T_i, scaled down, just to stay physically motivated
    # and numerically nonzero (some of your physics functions may not
    # like exactly n=0, e.g. sqrt/log terms in pressure or viscosity).
    n_small = 1e-6 * n_i

    bc_gas_A = n(x_A, t) ~ n_small

    return [ic_T, ic_n, bc_surface, bc_heat_B, bc_gas_B, bc_gas_A]

end
=#