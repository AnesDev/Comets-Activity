using ModelingToolkit
include("algorithm.jl")

function heliocentric_distance(t, a, e, M0)
    P = sqrt(a^3) * 365.25 * 24 * 60 * 60
    M = M0 + 2π * t / P
    E, _ = newton_raphson(M, e)
    return a * (1 - e * cos(E))
end

function initial_mean_anomaly(r0, a, e; inbound = true)
    E0 = acos((1 - r0 / a) / e)
    inbound && (E0 = -E0)
    return E0 - e * sin(E0)
end

const ORBIT = Ref((a = 1.0, e = 0.0, M0 = 0.0))

r_H_of_t(t::Number) = heliocentric_distance(t, ORBIT[].a, ORBIT[].e, ORBIT[].M0)
r_H_of_t(t::AbstractArray) = r_H_of_t.(t)
@register_symbolic r_H_of_t(t)