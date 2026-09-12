using ModelingToolkit
using DomainSets
using NeuralPDE

function build_problem(params, r_H, T_i_override, t_final)

    x_min = 0.0
    x_max = params.mantle_thickness
    t_min = 0.0
    t_max = t_final

    domains = [
        x ∈ Interval(0.0, x_max),
        t ∈ Interval(0.0, t_max)
    ]

    bcs, N_scale, T_i, n_i_scaled, n_A_scaled, solar_flux_scale, heat_flux_scale, gas_flux_scale =
        build_boundary_conditions(
            params,
            (x_A = x_min, x_B = x_max, t0 = t_min, t1 = t_max),
            r_H;
            T_i_override = T_i_override
        )

    eq_heat, eq_gas, heat_eq_scale, gas_rate_scale = build_equations(params, N_scale, T_i)
    equations = [eq_heat, eq_gas]

    println("eq types: ", typeof.(equations))
    println("bc types: ", typeof.(bcs))

    @named system = PDESystem(
        equations, bcs, domains, [x, t], [T(x, t), n(x, t)]
    )

    return system, N_scale, T_i, n_i_scaled, n_A_scaled, x_min, x_max, t_min, t_max,
           heat_eq_scale, gas_rate_scale, solar_flux_scale, heat_flux_scale, gas_flux_scale
end