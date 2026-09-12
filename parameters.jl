module Parameters

export Params, default_parameters

Base.@kwdef struct Params

    cond_A::Float64 = 1e-2
    cond_B::Float64 = 1e-3

    # conductivity
    sigma_SB::Float64 = 5.67e-8

    # Universal constants
    k_B::Float64 = 1.380649e-23
    R::Float64   = 8.314462618

    # Comet properties
    density::Float64 = 533.0
    porosity::Float64 = 0.80

    
    mantle_thickness::Float64 = 1.0
    t_final::Float64 = 1e4

    # Geometry
    pore_radius::Float64 = 1e-4
    CO_collision_diameter::Float64 = 3.76e-10

    # CO
    CO_molecular_mass::Float64 =
        28.01e-3 / 6.02214076e23


    F_sun::Float64  = 1360.8   # solar constant, W/m^2
    A_bond::Float64 = 0.04     # Bond albedo
    eps_IR::Float64 = 0.95     # emissivity
    T_OC::Float64   = 10.0     # Oort cloud / background temp, K

    activation_energy_CO::Float64 = 7450.0

    A_CO::Float64 = 8e9

    # Latent heat
    latent_heat_CO::Float64 = 7450.0 / 0.02801

end

default_parameters() = Params()

end