using Plots

"""
Plots the predicted temperature field.

Arguments
---------
x : Spatial coordinates.
t : Time coordinates.
T_pred : Predicted temperature matrix.
"""
function plot_temperature(x, t, T_pred)

    heatmap(
        x,
        t,
        T_pred,
        xlabel = "Depth (m)",
        ylabel = "Time (s)",
        title = "Temperature",
        colorbar_title = "K"
    )

end


"""
Plots the predicted CO molecular density.

Arguments
---------
x : Spatial coordinates.
t : Time coordinates.
n_pred : Predicted molecular density matrix.
"""
function plot_density(x, t, log10_n_pred)
    heatmap(
        x,
        t,
        log10_n_pred,
        xlabel = "Depth (m)",
        ylabel = "Time (s)",
        title = "CO Number Density",
        colorbar_title = "log₁₀(n) [m⁻³]"
    )
end

"""
Plots the optimization loss history.
"""
function plot_loss(losses)

    plot(
        losses,
        xlabel = "Iteration",
        ylabel = "Loss",
        yscale = :log10,
        title = "Training Loss",
        legend = false
    )

end