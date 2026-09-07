include("sublimation_pressure.jl")

"""
Upper-limit CO sublimation mass flux, Eq. (4):
    f_CO = P(T) * sqrt(m / (2 pi k_B T))
Units: kg m^-2 s^-1.
"""
function sublimation_flux(T, params)
    P = sublimation_pressure(T, params)
    return P * sqrt(params.CO_molecular_mass / (2 * pi * params.k_B * T))
end