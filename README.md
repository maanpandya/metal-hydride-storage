# Numerical Simulation of Hydrogen Storage in Metal-Hydride Tanks

![Julia](https://img.shields.io/badge/Julia-1.9%2B-blueviolet)
![License](https://img.shields.io/badge/License-MIT-green)
![Status](https://img.shields.io/badge/Status-Active-blue)

This repository contains the collaborative work for the TAship focused on the mathematical modeling and numerical simulation of hydrogen storage in metal-oxide frameworks. Our goal is to develop and analyze models for the laminar flow, absorption, and desorption of hydrogen gas under iso-thermal conditions.

## 📝 Project Overview

This project aims to simulate the complex physics involved in hydrogen storage tanks that use metal-oxide frameworks as a porous binding agent. The work is divided into several key areas:
1.  **Fluid Dynamics:** Modeling the laminar flow of hydrogen gas through the porous medium using formulations like Stokes and Navier-Stokes.
2.  **Reaction Kinetics:** Simulating the absorption (binding) and desorption (release) of hydrogen from the metal-oxide material.
3.  **Numerical Methods:** Implementing these models using the finite element method (FEM) and advanced time-integration techniques.

The primary software stack for this project is built on the Julia programming language, leveraging its high-performance scientific computing ecosystem.

## 🛠️ Technology Stack

*   **Core Language:** [Julia](https://julialang.org/)
*   **Finite Element Method (FEM):** [Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl/stable/)
*   **Meshing:** [GMSH.jl](https://github.com/JuliaFEM/Gmsh.jl) for geometry definition and mesh generation.
*   **Time Integration (ODEs/DAEs):** [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/)
*   **Numerical Linear Algebra:** [LinearSolve.jl](https://docs.sciml.ai/LinearSolve/stable/), [KrylovKit.jl](https://jutho.github.io/KrylovKit.jl/stable/), [Preconditioners.jl](https://julialinearalgebra.github.io/Preconditioners.jl/stable/)
*   **Visualization:** [Plots.jl](https://docs.juliaplots.org/stable/), [Paraview](https://www.paraview.org/) (for `.vtu` files)

## 📂 Repository Structure

To keep our work organized, we will follow this directory structure:

TBD


## 🚀 Getting Started

To set up your local environment and run the code, follow these steps:

1.  **Clone the repository:**
    ```bash
    git clone www.github.com/maanpandya/metal-hydride-storage
    cd <repository-directory>
    ```

2.  **Install Julia:**
    Make sure you have Julia installed (version 1.9 or later is recommended). You can download it from the [official Julia website](https://julialang.org/downloads/).

3.  **Activate the Julia Environment:**
    Open a Julia REPL in the project's root directory and activate the shared environment. This will ensure everyone uses the same package versions.
    ```julia
    julia> ]  # Press ']' to enter Pkg mode
    pkg> activate .
    pkg> instantiate
    ```
    The `instantiate` command will download and install all the necessary packages listed in `Project.toml` and `Manifest.toml`.

## 🧑‍💻 Team Roles and Responsibilities

Our work is divided based on a "Divide and Conquer" approach. Each TA has a primary area of responsibility.

### **TA-1: Swayam Kuckreja (Meshing & Visualization)**

Swayam is responsible for all pre-processing (geometry and meshing) and post-processing (visualization) tasks. This provides a stable foundation for the simulation team.
*   **Geometry & Meshing:** Manage geometry definition and mesh generation using `GMSH.jl`, creating 2D and 3D meshes for the tank and nozzle.
*   **Boundary Labeling:** Ensure all mesh patches (inlets, walls, outlets) are correctly labeled for assigning boundary conditions.
*   **Mesh Integration:** Verify that meshes can be correctly read by `Ferrite.jl`.
*   **Visualization:** Support the team in visualizing meshes and simulation results using `Paraview` and other tools.

### **TA-2: Maan Pandya (Time Integration & Transient Modeling)**

Maan is responsible for setting up and solving time-dependent problems, bridging the 0D models with the full FEM simulations.
*   **0D Modeling:** Implement and solve 0D reactor models for pressure, temperature, and species density using `DifferentialEquations.jl`.
*   **Transient FEM:** Run and analyze transient Stokes and Navier-Stokes models from `Ferrite.jl` tutorials.
*   **Solver Profiling:** Analyze the performance of different time-integration methods and solver settings.
*   **Advanced Solvers:** Explore and implement advanced iterative linear solvers (`GMRES`, `ILU`) to accelerate implicit time-stepping.

### **TA-3: Rimaz Khan (FEM Prototyping & Physics Integration)**

Rimaz is responsible for the rapid prototyping of new physics within the `Ferrite.jl` framework, using simpler, internally generated meshes.
*   **Physics Formulation:** Extend the Stokes/Navier-Stokes models to include coupled transport equations for H2 gas and metal-oxide densities.
*   **Jacobian and Residual Extension:** Modify the core `Ferrite.jl` assembly loops to expand the system from a 2x2 (pressure-velocity) to a 4x4 block system.
*   **Weak & Strong Coupling:** Implement and test both weakly coupled (post-processing) and strongly coupled (fully implicit) physics integrations.

### **TA-4: Simranjeet Singh (FEM Implementation on Complex Geometries)**

Simranjeet is responsible for implementing the advanced physics models developed by Rimaz on the complex, realistic geometries provided by Swayam.
*   **Model Porting:** Adapt the prototype models from rectangular channels to the full tank-with-nozzle geometries.
*   **Mesh Integration:** Work closely with Swayam to use the `GMSH.jl`-generated meshes in `Ferrite.jl` simulations.
*   **Boundary Conditions:** Implement the correct boundary conditions on the complex geometries.
*   **Scaling Up:** Run the final, fully-featured simulations on fine meshes to produce high-fidelity results.

## 🤝 Contribution Guidelines

To maintain a clean and collaborative workflow, please follow these guidelines:

1.  **Branching:** Create a new branch for every new feature or task. Name it descriptively, e.g., `feature/maan-0d-model` or `fix/swayam-mesh-labels`.
2.  **Commits:** Write clear and concise commit messages.
3.  **Pull Requests (PRs):** When your feature is complete, open a Pull Request to merge your branch into `main`. Briefly describe the changes in the PR.
4.  **Code Style:** Please format your Julia code using [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) to maintain a consistent style across the project.

## ⚖️ License

This project is licensed under the MIT License. See the `LICENSE` file for details.
