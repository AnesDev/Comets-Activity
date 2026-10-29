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

include("Orbit/comet_orbit.jl")
include("equations.jl")
include("boundary_conditions.jl")
include("problem.jl")
include("discretisation.jl")
include("plots.jl")
include("train.jl")

mkpath("checkpoints")
mkpath("results")

const A_ORB, E_ORB = 51.0, 0.98

cases = [
    (r_H0 = 100.0, t_final = 2.5e9,  label = "100au"),   # ends at ≈ 80 AU
    (r_H0 = 10.0,  t_final = 4.4e7,  label = "10au"),    # ends at ≈ 6 AU
    (r_H0 = 4.0,   t_final = 1.2e7,  label = "4au"),     # ends at ≈ 2.4 AU
]

for case in cases
    println("========== Running r_H0 = $(case.r_H0) AU ==========")

    ORBIT[] = (a = A_ORB, e = E_ORB,
               M0 = initial_mean_anomaly(case.r_H0, A_ORB, E_ORB; inbound = true))

    params = Parameters.Params(t_final = case.t_final)

    system, N_scale, T_i, n_i_scaled, n_A_scaled, x_min, x_max, t_min, t_max,
    heat_eq_scale, gas_rate_scale, solar_flux_scale, heat_flux_scale, gas_flux_scale =
        build_problem(params, case.r_H0, nothing, case.t_final)

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