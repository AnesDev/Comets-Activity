function knudsen_diffusion(T, params)

    return thermal_speed(T, params) *
           params.pore_radius / 3

end