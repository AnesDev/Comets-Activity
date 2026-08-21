using NeuralPDE

include("model.jl")

model = build_model(params)   # params must already exist — see main.jl order below

strategy = QuasiRandomTraining(2000)

discretization = PhysicsInformedNN(
    model,
    strategy
)