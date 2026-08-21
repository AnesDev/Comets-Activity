"""
Ideal gas pressure.
"""
function pressure(n, T, params)
    return n * params.k_B * T
end