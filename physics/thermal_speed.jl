"""
Thermal speed of CO molecules.

Arguments
---------
T : Temperature [K]

Returns
-------
Thermal speed [m/s]
"""
function thermal_speed(T, params)

    return sqrt(
        8 * params.k_B * T /
        (π * params.CO_molecular_mass)
    )

end