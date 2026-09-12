using Lux

function build_model(params, x_A, x_B, T_i, n_i_scaled, n_A_scaled)

    x_scale = params.mantle_thickness
    t_scale = params.t_final
    T_floor = params.T_OC

    normalize_inputs(xt) = vcat(
        xt[1:1, :] ./ x_scale,
        xt[2:2, :] ./ t_scale
    )

    core = Chain(
        normalize_inputs,
        Dense(2, 64, tanh),
        Dense(64, 64, tanh),
        Dense(64, 64, tanh),
        Dense(64, 2),
    )

    function scale_outputs(core_out, xt)
        T_raw = core_out[1:1, :]
        n_raw = core_out[2:2, :]

        x = xt[1:1, :]
        t = xt[2:2, :]
        envelope_t = t ./ t_scale

        T_scaled = (1 .- envelope_t) .* T_i .+ envelope_t .* (T_floor .+ abs.(T_raw))

        n_free = (1 .- envelope_t) .* n_i_scaled .+ envelope_t .* softplus.(n_raw)

        envelope_x = (x .- x_A) ./ (x_B .- x_A)

        n_scaled = n_A_scaled .+ envelope_x .* (n_free .- n_A_scaled)

        return vcat(T_scaled, n_scaled)
    end

    return Chain(SkipConnection(core, scale_outputs))

end