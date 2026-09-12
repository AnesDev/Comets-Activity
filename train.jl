using NeuralPDE
using Optimization
using OptimizationOptimisers
using OptimizationOptimJL
using Serialization
using Dates

function run_training(system, discretization, N_scale, x_min, x_max, t_min, t_max,
                       checkpoint_file, label,
                       heat_eq_scale, gas_rate_scale, solar_flux_scale, heat_flux_scale, gas_flux_scale;
                       checkpoint_every=100)

    losses = Float64[]
    start_time = time()
    phi = discretization.phi

    log_file = replace(checkpoint_file, ".jls" => "_log.txt")
    io = open(log_file, "w")

    function logprintln(args...)
        println(args...)
        println(io, args...)
        flush(io)
    end

    watch_points = [
        ("IC  @ x=A, t=0", x_min, t_min),
        ("IC  @ x=B, t=0", x_max, t_min),
        ("BC  @ x=B, t=t_max", x_max, t_max),
    ]

    function print_watch_points(u)
        for (label_, xv, tv) in watch_points
            out = phi([xv, tv], u)
            T_val = out[1]
            n_val = N_scale * out[2]
            logprintln("    $label_ -> T=$(round(T_val, sigdigits=5)) K, n=$(round(n_val, sigdigits=5)) m^-3")
        end
    end

    sym_prob = symbolic_discretize(system, discretization)

    scale_info = [
        ("PDE[1] heat eq", heat_eq_scale,   "W/m^3"),
        ("PDE[2] gas eq",  gas_rate_scale,  "m^-3/s"),
        ("BC[1] surface",  solar_flux_scale,"W/m^2"),
        ("BC[2] ice heat", heat_flux_scale, "W/m^2"),
        ("BC[3] ice gas",  gas_flux_scale,  "m^-2/s"),
    ]

    function print_loss_breakdown(u)
        logprintln("========== LOSS BREAKDOWN ==========")
        pde_losses = map(f -> f(u), sym_prob.loss_functions.pde_loss_functions)
        bc_losses  = map(f -> f(u), sym_prob.loss_functions.bc_loss_functions)
        all_losses = vcat(pde_losses, bc_losses)

        for (i, (name, scale, unit)) in enumerate(scale_info)
            loss_val = all_losses[i]
            physical = sqrt(loss_val) * scale
            logprintln("$name : scaled=$(loss_val)   physical≈$(physical) $unit")
        end

        logprintln("Total PDE = ", sum(pde_losses))
        logprintln("Total BC  = ", sum(bc_losses))
        logprintln("Grand Total = ", sum(pde_losses) + sum(bc_losses))
        logprintln("===============================")
    end

    callback = function (state, loss)
        push!(losses, loss)
        elapsed = time() - start_time
        logprintln("[$(length(losses))] loss = $(round(loss, sigdigits=6))  ($(round(elapsed, digits=1))s)")
        if length(losses) % 20 == 1
            print_loss_breakdown(state.u)
        end
        print_watch_points(state.u)

        if length(losses) % checkpoint_every == 0
            iter_checkpoint_file = replace(checkpoint_file, ".jls" => "_iter$(length(losses)).jls")
            Serialization.serialize(iter_checkpoint_file, state.u)
            logprintln("checkpoint saved at iteration $(length(losses)) -> $iter_checkpoint_file")
            notify_checkpoint(label, length(losses))
        end

        if length(losses) % 100 == 0
            notify_progress(label, length(losses), loss, elapsed)
            notify_loss_plot(label, length(losses), loss, losses)
        end

        return false
    end

    prob = discretize(system, discretization)

    logprintln()
    logprintln("========== LOSS FUNCTIONS ==========")
    logprintln("PDE losses:")
    for (i, _) in enumerate(sym_prob.loss_functions.pde_loss_functions)
        logprintln("  PDE $i")
    end
    logprintln("BC losses:")
    logprintln("BC[1] = Surface Energy")
    logprintln("BC[2] = Ice Heat Flux")
    logprintln("BC[3] = Ice Gas Flux")
    logprintln("====================================")
    logprintln()

    result = Optimization.solve(prob, OptimizationOptimisers.Adam(1e-3), callback=callback, maxiters=600)
    Serialization.serialize(checkpoint_file, result.u)

    prob = remake(prob; u0 = result.u)
    try
        result = Optimization.solve(prob, LBFGS(), callback=callback, maxiters=400)
    catch err
        logprintln("LBFGS stopped early: $err")
        notify_error(label, err)
    end
    Serialization.serialize(checkpoint_file, result.u)

    notify_finish(label, losses[end], time() - start_time)

    close(io)

    return result, losses, phi
end