using NeuralPDE

include("model.jl")

model = build_model(params, x_min, x_max, T_i, n_i_scaled, n_A_scaled)

strategy = QuasiRandomTraining(2000)


discretization = PhysicsInformedNN(
    model,
    strategy
    )