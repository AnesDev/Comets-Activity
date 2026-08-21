const ALPHA = 1.0

"""
Dynamic viscosity of CO gas:
    μ = (α / π^(3/2)) * sqrt(k_B m T) / σ²

Arguments
---------
T : Temperature [K]

Returns
-------
Dynamic viscosity [Pa·s]
"""
function viscosity(T, params)
    return (ALPHA / π^(3/2)) *
           sqrt(params.k_B * params.CO_molecular_mass * T) /
           params.CO_collision_diameter^2
end