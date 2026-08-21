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

include("problem.jl")        # <-- moved up: this creates `params`
include("discretisation.jl") # <-- now calls build_model(params) safely, params exists

include("plots.jl")
include("discord.jl")

include("train.jl")