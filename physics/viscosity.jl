const ALPHA = 1.0

"""
Dynamic viscosity of CO gas
"""
function viscosity(T, params)
    return (ALPHA / π^(3/2)) *
           sqrt(params.k_B * params.CO_molecular_mass * T) /
           params.CO_collision_diameter^2
end