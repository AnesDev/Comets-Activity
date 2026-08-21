"""
Permeability of the porous medium.

Returns
-------
Permeability [m²]
"""
function permeability(params)

    return params.pore_radius^2 / 32

end