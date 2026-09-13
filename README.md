# Comets-Activity

A Physics-Informed Neural Network (PINN) framework for modeling comet activity and gas sublimation dynamics across different heliocentric distances.

## Overview

This internship project develops a neural network-based solver for simulating the thermal and gas diffusion processes occurring in cometary mantles. The framework solves coupled partial differential equations (PDEs) governing heat conduction and gas transport in cometary nuclei, with applications across multiple heliocentric distances (4 AU, 10 AU, and 100 AU from the Sun).

### Key Features

- **Physics-Informed Neural Networks**: Integrates domain knowledge through differential equations directly into the learning process
- **Multi-scale Simulations**: Handles different heliocentric distances with appropriate physical parameters
- **Coupled PDE System**: Simultaneously solves heat transport and gas diffusion equations
- **Automatic Differentiation**: Leverages Julia's Enzyme and ForwardDiff for computing derivatives
- **Visualization & Tracking**: Includes real-time training monitoring with Discord notifications

## Repository Structure

```
Comets-Activity/
├── main.jl                      # Entry point - runs training pipeline for all cases
├── equations.jl                 # PDE system definition using ModelingToolkit
├── boundary_conditions.jl       # Boundary conditions for heat and gas equations
├── model.jl                     # Neural network architecture (Lux-based)
├── train.jl                     # Training loop and optimization
├── discretisation.jl            # Domain discretization strategy
├── parameters.jl                # Physical constants and model parameters
├── plots.jl                     # Visualization functions for results
├── discord.jl                   # Discord notifications for training events
├── physics/                     # Physics submodule
│   ├── pressure.jl
│   ├── heat_capacity.jl
│   ├── conductivity.jl
│   ├── thermal_speed.jl
│   ├── viscosity.jl
│   ├── permeability.jl
│   ├── knudsen_diffusion.jl
│   ├── sublimation_pressure.jl
│   └── sublimation_flux.jl
├── checkpoints/                 # Training checkpoints (auto-created)
├── results/                     # Output plots and results (auto-created)
├── Project.toml                 # Julia package dependencies
└── Manifest.toml                # Locked dependency versions
```

## Physical Model

### Governing Equations

The framework solves a coupled system of PDEs:

#### Heat Equation
```
ρ c_p ∂T/∂t = ∂/∂x(κ(T) ∂T/∂x)
```

Where:
- **ρ**: Bulk density of cometary material
- **c_p**: Heat capacity (temperature-dependent)
- **κ(T)**: Thermal conductivity (temperature-dependent)
- **T**: Temperature

#### Gas Transport Equation
```
φ ∂n/∂t = ∂/∂x(D_K(T) ∂n/∂x + (k/η(T)) n k_B (∂n/∂x T + n ∂T/∂x))
```

Where:
- **φ**: Porosity of the mantle
- **n**: Gas density
- **D_K**: Knudsen diffusion coefficient
- **η(T)**: Gas viscosity (temperature-dependent)
- **k_B**: Boltzmann constant

### Physical Processes

The model incorporates realistic physics:

- **Sublimation**: CO ice sublimation driven by solar heating
- **Heat Conduction**: Temperature-dependent thermal conductivity through the mantle
- **Gas Diffusion**: Combined Knudsen diffusion and pressure-driven flow
- **Temperature-Dependent Properties**: Heat capacity, viscosity, and conductivity vary with T
- **Molecular Transport**: Molecular mass, thermal speeds, and permeability effects

## Simulation Cases

The framework runs three benchmark cases at different heliocentric distances:

| Case | Heliocentric Distance | Initial Temperature | Simulation Time |
|------|-----------------------|---------------------|-----------------|
| 100 AU | 100 AU | 40 K | 10,000 seconds |
| 10 AU | 10 AU | Computed | 150 seconds |
| 4 AU | 4 AU | Computed | 100 seconds |

Each case produces:
- **Temperature Profile**: T(x,t) solution field
- **Density Profile**: n(x,t) solution field
- **Loss Curve**: Training loss evolution

## Neural Network Architecture

The physics-informed neural network uses:

- **Input Layer**: Normalized spatial (x) and temporal (t) coordinates
- **Hidden Layers**: 3 × 64-neuron layers with tanh activations
- **Output Layer**: 2 outputs (temperature T and gas density n)
- **Output Scaling**: Physics-consistent envelope functions that enforce:
  - Initial and boundary conditions
  - Physical monotonicity constraints
  - Proper scaling with domain extent

### Key Features

- Input normalization by domain scales
- Output scaling via learnable envelopes
- Skip connections for improved training
- Softplus activation for ensuring positivity of density

## Dependencies

Core scientific libraries:

- **[NeuralPDE.jl](https://github.com/SciML/NeuralPDE.jl)**: Physics-informed neural networks
- **[ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl)**: Symbolic PDE definition
- **[Lux.jl](https://github.com/LuxDL/Lux.jl)**: Neural network framework
- **[Optimization.jl](https://github.com/SciML/Optimization.jl)**: Optimization algorithms
- **[Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl)**: Automatic differentiation

Utilities:

- **CairoMakie.jl**: Visualization backend
- **ComponentArrays.jl**: Flexible array structures
- **ForwardDiff.jl**: Forward-mode differentiation
- **JSON3.jl**: Data serialization
- **HTTP.jl**: Network communication for Discord notifications

See `Project.toml` for complete dependency list.

## Installation & Setup

### Prerequisites

- Julia ≥ 1.9
- Git

### Steps

1. **Clone the repository**:
   ```bash
   git clone https://github.com/AnesDev/Comets-Activity.git
   cd Comets-Activity
   ```

2. **Activate the Julia environment**:
   ```bash
   julia --project
   ```

3. **Instantiate dependencies**:
   ```julia
   julia> ]instantiate
   ```

4. **(Optional) Configure Discord notifications**:
   Edit `discord.jl` to set your webhook URL:
   ```julia
   WEBHOOK_URL = "your-discord-webhook-url"
   ```

## Running Simulations

### Standard Run (All Cases)

```bash
julia main.jl
```

This runs training for all three heliocentric cases (100 AU, 10 AU, 4 AU) sequentially.

### Output Files

After completion, check:

- **`results/temperature_*.png`**: Temperature profiles across space-time
- **`results/density_*.png`**: Log₁₀ gas density profiles
- **`results/loss_*.png`**: Training loss curves
- **`checkpoints/checkpoint_*.jls`**: Saved model parameters (every 100 epochs)

## File Descriptions

### Core Components

| File | Purpose |
|------|---------|
| `main.jl` | Orchestrates the complete pipeline for all simulation cases |
| `equations.jl` | Defines the coupled PDE system using symbolic modeling |
| `model.jl` | Implements the neural network architecture |
| `train.jl` | Training loop, loss computation, and optimization |
| `boundary_conditions.jl` | Specifies BCs for heat and gas equations |
| `problem.jl` | Assembles the complete problem for NeuralPDE |
| `discretisation.jl` | Configures spatial/temporal mesh strategy |
| `parameters.jl` | Physical constants and customizable parameters |

### Utilities

| File | Purpose |
|------|---------|
| `plots.jl` | Functions for visualizing T, n, and loss |
| `discord.jl` | Sends training status/results to Discord |
| `physics/*.jl` | Submodules for physical property calculations |

## Training Details

### Optimization Strategy

- **Algorithm**: Optimized using `Optimization.jl` (default: Adam + L-BFGS)
- **Loss Function**: Physics residual minimization (PDE satisfaction)
- **Checkpointing**: Models saved every 100 iterations
- **Automatic Differentiation**: Enzyme-based automatic differentiation

### Scaling & Normalization

- **Length Scale**: Mantle thickness
- **Time Scale**: Simulation final time
- **Temperature Scale**: Heat flux driven by sublimation
- **Density Scale**: Sublimation-driven gas flux

This ensures well-conditioned optimization problems.

## Expected Results

### Temperature Solutions

- Rapid heating in early times near surface
- Thermal diffusion propagating into the mantle
- Equilibration over the simulation duration
- Stronger heating at smaller heliocentric distances (closer to Sun)

### Density Solutions

- Initial gas accumulation from sublimation
- Diffusive spreading into porous mantle
- Competition between Knudsen diffusion and pressure-driven flow
- More pronounced gradients at smaller distances

### Convergence

Typical loss reduction: 10⁻¹ → 10⁻⁵ over training (case-dependent)

## Project Context

This is an internship project focused on developing and validating Physics-Informed Neural Networks for cometary science applications. The work bridges:

- **Applied Mathematics**: PINN methodology and optimization
- **Physics**: Cometary physics, sublimation, and transport
- **Scientific Computing**: Large-scale neural network training
- **Software Engineering**: Modular Julia code design

