using NeuralPDE
using Optimization
using OptimizationOptimisers
using OptimizationOptimJL
using Serialization
using Dates

const CHECKPOINT_FILE = "checkpoint.jls"
const CHECKPOINT_EVERY = 100

notify_start("Comet Activity", strategy)

losses = Float64[]
start_time = time()

phi = discretization.phi

# Points worth watching: surface at t=0, ice front at t=0, and ice
# front at t=t_max (to see how the BC-driven behavior evolves over
# the whole training window, not just at t=0).
watch_points = [
    ("IC  @ x=A, t=0", x_min, t_min),
    ("IC  @ x=B, t=0", x_max, t_min),
    ("BC  @ x=B, t=t_max", x_max, t_max),
]

function print_watch_points(u)
    for (label, xv, tv) in watch_points
        out = phi([xv, tv], u)
        T_val = out[1]
        n_val = N_scale * out[2]
        println("    $label -> T=$(round(T_val, sigdigits=5)) K, n=$(round(n_val, sigdigits=5)) m^-3")
    end
end

function print_loss_breakdown(u)

    println("========== LOSS BREAKDOWN ==========")

    pde_losses = map(f -> f(u), sym_prob.loss_functions.pde_loss_functions)
    bc_losses  = map(f -> f(u), sym_prob.loss_functions.bc_loss_functions)

    for (i,l) in enumerate(pde_losses)
        println("PDE[$i] = $(l)")
    end

    for (i,l) in enumerate(bc_losses)
        println("BC[$i] = $(l)")
    end

    println("Total PDE = ", sum(pde_losses))
    println("Total BC  = ", sum(bc_losses))
    println("Grand Total = ", sum(pde_losses) + sum(bc_losses))

    println("===============================")
end

callback = function (state, loss)

    push!(losses, loss)

    elapsed = time() - start_time

    println("[$(length(losses))] loss = $(round(loss, sigdigits=6))  ($(round(elapsed, digits=1))s)")
    if length(losses) % 20 == 1
        print_loss_breakdown(state.u)
    end
    print_watch_points(state.u)

    if length(losses) % CHECKPOINT_EVERY == 0
        Serialization.serialize(CHECKPOINT_FILE, state.u)
        println("checkpoint saved at iteration $(length(losses))")
        notify_checkpoint(length(losses))
    end

    if length(losses) % 100 == 0

        notify_progress(
            length(losses),
            loss,
            elapsed
        )

        notify_loss_plot(
            length(losses),
            loss,
            losses
        )

    end

    return false

end


sym_prob = symbolic_discretize(system, discretization)

prob = discretize(system, discretization)

println()
println("========== LOSS FUNCTIONS ==========")

println("PDE losses:")
for (i, _) in enumerate(sym_prob.loss_functions.pde_loss_functions)
    println("  PDE $i")
end

println("BC losses:")
println("BC[1] = IC Temperature")
println("BC[2] = IC Density")
println("BC[3] = Surface Energy")
println("BC[4] = Ice Heat Flux")
println("BC[5] = Ice Gas Flux")
println("BC[6] = Surface Vacuum")

println("====================================")
println()


result = Optimization.solve(
    prob,
    OptimizationOptimisers.Adam(1e-3),
    callback = callback,
    maxiters = 600
)

Serialization.serialize(CHECKPOINT_FILE, result.u)

prob = remake(prob; u0 = result.u)

result = Optimization.solve(
    prob,
    LBFGS(),
    callback = callback,
    maxiters = 400
)

Serialization.serialize(CHECKPOINT_FILE, result.u)

notify_finish(
    losses[end],
    time() - start_time
)

xs = range(x_min, x_max, length = 100)
ts = range(t_min, t_max, length = 100)

T_pred = [phi([x, t], result.u)[1] for t in ts, x in xs]
n_pred = [N_scale * phi([x, t], result.u)[2] for t in ts, x in xs]

p_T = plot_temperature(xs, ts, T_pred)
p_n = plot_density(xs, ts, log10.(n_pred))
p_loss = plot_loss(losses)

savefig(p_T, "temperature.png")
savefig(p_n, "density.png")
savefig(p_loss, "loss.png")

notify_results(
    solution = "temperature.png",
    density = "density.png",
    loss = "loss.png"
)