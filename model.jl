using Lux

function build_model(params)

    x_scale = params.mantle_thickness
    t_scale = params.t_final

    normalize_inputs(xt) = vcat(
        xt[1:1, :] ./ x_scale,
        xt[2:2, :] ./ t_scale
    )

    function scale_outputs(x)
        T_raw = x[1:1, :]
        n_raw = x[2:2, :]
        T_scaled = params.T_OC .+ abs.(T_raw)   # physical floor at T_OC, constant-strength gradient
        n_scaled = softplus.(n_raw)
        return vcat(T_scaled, n_scaled)
    end

    return Chain(
        normalize_inputs,
        Dense(2, 64, tanh),
        Dense(64, 64, tanh),
        Dense(64, 64, tanh),
        Dense(64, 2),
        scale_outputs
    )

end