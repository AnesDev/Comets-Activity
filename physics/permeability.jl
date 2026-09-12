"""
Permeability of the porous medium.
"""
function permeability(params)

    return params.pore_radius^2 / 32

end