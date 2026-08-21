"""
Thermal conductivity of the porous comet material, Eq. (2):
    k(T) = A + B*(4*eps_IR*sigma_SB*T^3)

Arguments
---------
T      : Temperature [K]
params : Params struct (uses cond_A, cond_B, eps_IR, sigma_SB)

Returns
-------
Thermal conductivity [W/(m·K)]
"""
function conductivity(T, params)
    return params.cond_A +
           params.cond_B *
           (4 * params.eps_IR * params.sigma_SB * T^3)
end