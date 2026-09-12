using NeuralPDE
include("model.jl")

function build_discretization(params, x_min, x_max, T_i, n_i_scaled, n_A_scaled)
    model = build_model(params, x_min, x_max, T_i, n_i_scaled, n_A_scaled)
    strategy = QuasiRandomTraining(2000)
    return PhysicsInformedNN(model, strategy)
end