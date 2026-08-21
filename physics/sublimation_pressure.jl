"""
Equilibrium sublimation pressure of CO ice, Eq. (5):
    P(T) = A_CO * exp(-E_CO / (R T))
Fitted values from Table 1: A_CO = 8e9 N/m^2, activation energy = 7450 J/mol.
"""
function sublimation_pressure(T, params)
    return params.A_CO * exp(-params.activation_energy_CO / (params.R * T))
end