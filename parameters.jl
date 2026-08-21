module Parameters

export Params, default_parameters

Base.@kwdef struct Params

    cond_A::Float64 = 1e-2   # W/(m·K) — sweep per your supervisor's suggestion
    cond_B::Float64 = 1e-3   # dimensionless attenuation factor

    # conductivity
    sigma_SB::Float64 = 5.67e-8

    # Universal constants
    k_B::Float64 = 1.380649e-23
    R::Float64   = 8.314462618

    # Comet properties
    density::Float64 = 533.0
    porosity::Float64 = 0.80

    
    mantle_thickness::Float64 = 1.0    # meters — this is x_B in boundary_conditions.jl
    t_final::Float64 = 1e4             # seconds — placeholder; see note below

    # Geometry
    pore_radius::Float64 = 1e-4
    CO_collision_diameter::Float64 = 3.76e-10

    # CO
    CO_molecular_mass::Float64 =
        28.01e-3 / 6.02214076e23

    # --- Needed for boundary_conditions.jl ---

    # Orbital / radiative (Eq. 1, Table 1)
    F_sun::Float64  = 1360.8   # solar constant, W/m^2
    A_bond::Float64 = 0.04     # Bond albedo
    eps_IR::Float64 = 0.95     # emissivity
    T_OC::Float64   = 10.0     # Oort cloud / CMB background temp, K

    # --- CO sublimation: TWO separate quantities, NOT the same ---
    # Table 1 reuses one symbol (E_CO) for both, but the units differ
    # and the paper's own Appendix confirms they're different numbers
    # ("L = 2.66e5 J/kg is the latent heat of CO from Table 1").

    # Arrhenius activation energy, Eq. (5), used inside P(T). Units: J/mol.
    activation_energy_CO::Float64 = 7450.0

    # Arrhenius prefactor, Eq. (5). Units: N/m^2.
    A_CO::Float64 = 8e9

    # Latent heat, Eq. (3), used in the heat BC at x=B. Units: J/kg.
    # = activation_energy_CO / molar_mass_CO ≈ 2.66e5 J/kg (matches
    # the Appendix's stated value almost exactly).
    latent_heat_CO::Float64 = 7450.0 / 0.02801

end

default_parameters() = Params()

end