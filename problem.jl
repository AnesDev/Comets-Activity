using ModelingToolkit
using DomainSets
using NeuralPDE

params = Parameters.default_parameters()

########################################
# Domain
########################################

x_min = 0.0
x_max = params.mantle_thickness

t_min = 0.0
t_max = params.t_final

domains = [
    x ∈ Interval(0.0, x_max),
    t ∈ Interval(0.0, t_max)
]

########################################
# Boundary / Initial Conditions (also gives us N_scale)
########################################
bcs, N_scale = build_boundary_conditions(
    params,
    (x_A = x_min, x_B = x_max, t0 = t_min, t1 = t_max),
    100.0      # heliocentric distance (AU)
)

########################################
# Equations (needs N_scale from above)
########################################
eq_heat, eq_gas = build_equations(params, N_scale)
equations = [eq_heat, eq_gas]

########################################
# PDE System
########################################

println("eq types: ", typeof.(equations))
println("bc types: ", typeof.(bcs))

@named system = PDESystem(
    equations,
    bcs,
    domains,
    [x, t],
    [T(x, t), n(x, t)]
)