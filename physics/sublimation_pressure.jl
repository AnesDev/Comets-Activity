"""
Equilibrium sublimation pressure of CO ice
    P(T) = A_CO * exp(-E_CO / (R T))
"""
function sublimation_pressure(T, params)
    return params.A_CO * exp(-params.activation_energy_CO / (params.R * T))
end