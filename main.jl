using Random

include("parameters.jl")
using .Parameters

include("physics/pressure.jl")
include("physics/heat_capacity.jl")
include("physics/conductivity.jl")
include("physics/thermal_speed.jl")
include("physics/viscosity.jl")
include("physics/permeability.jl")
include("physics/knudsen_diffusion.jl")
include("physics/sublimation_pressure.jl")
include("physics/sublimation_flux.jl")

include("equations.jl")
include("boundary_conditions.jl")
include("problem.jl")
include("discretisation.jl")
include("plots.jl")
include("train.jl")

mkpath("checkpoints")
mkpath("results")

cases = [
    (r_H = 100.0, T_i = 40.0, t_final = 1e4,   label = "100au"),
    (r_H = 10.0,  T_i = nothing, t_final = 150.0, label = "10au"),
    (r_H = 4.0,   T_i = nothing, t_final = 100.0, label = "4au"),
]

for case in cases
    println("========== Running r_H = $(case.r_H) AU ==========")

    params = Parameters.default_parameters()

    system, N_scale, T_i, n_i_scaled, n_A_scaled, x_min, x_max, t_min, t_max,
    heat_eq_scale, gas_rate_scale, solar_flux_scale, heat_flux_scale, gas_flux_scale =
        build_problem(params, case.r_H, case.T_i, case.t_final)

    discretization = build_discretization(params, x_min, x_max, T_i, n_i_scaled, n_A_scaled)

    checkpoint_file = joinpath("checkpoints", "checkpoint_$(case.label).jls")
    result, losses, phi = run_training(
        system, discretization, N_scale, x_min, x_max, t_min, t_max,
        checkpoint_file, case.label,
        heat_eq_scale, gas_rate_scale, solar_flux_scale, heat_flux_scale, gas_flux_scale;
        checkpoint_every=100
    )

    xs = range(x_min, x_max, length=100)
    ts = range(t_min, t_max, length=100)
    T_pred = [phi([x, t], result.u)[1] for t in ts, x in xs]
    n_pred = [N_scale * phi([x, t], result.u)[2] for t in ts, x in xs]

    temp_file = joinpath("results", "temperature_$(case.label).png")
    dens_file = joinpath("results", "density_$(case.label).png")
    loss_file = joinpath("results", "loss_$(case.label).png")

    savefig(plot_temperature(xs, ts, T_pred), temp_file)
    savefig(plot_density(xs, ts, log10.(n_pred)), dens_file)
    savefig(plot_loss(losses), loss_file)

end