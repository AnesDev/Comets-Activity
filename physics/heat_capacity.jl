"""
Specific heat capacity of water ice [J/(kg·K)]
"""

function cp_ice(T)
    return 7.73e-3 * T *
           (1 - exp(-1.263e-3 * T^2)) *
           (
               1 +
               exp(-3 * sqrt(T)) * 8.47e-3 * T^6 +
                   2.0825e-7 * T^4 * exp(-4.97e-2 * T)
           )
end


"""
Specific heat capacity of cometary dust [J/(kg·K)]
"""
function cp_dust(T)

    return 2.03 * T *
           (1 - exp(-2.5e-4 * T^2)) *
           (
               1 +
               1.24e-9 * T^6 * exp(-0.853 * sqrt(T)) -
               3.69e6 * T^4 * exp(-4.76e5 * T)
           )

end

"""
Specific heat capacity of the comet mantle

The mantle is assumed to consist of a 1:1 mixture of
water ice and refractory dust.

Arguments
"""
function heat_capacity(T)

    return (cp_ice(T) + cp_dust(T)) / 2

end