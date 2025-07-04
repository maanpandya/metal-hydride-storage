# This version uses internal energy U as a state variable instead of temperature T,
# matching the physics of the original implementation.

using DifferentialEquations
using Plots

# Parameters struct
struct TankParameters
    P0::Float64     # Inlet total pressure [Pa]
    T0::Float64     # Inlet total temperature [K]
    V::Float64      # Tank volume [m³]
    gamma::Float64  # Specific heat ratio [-]
    R::Float64      # Specific gas constant [J/(kg·K)]
    Cd::Float64     # Discharge coefficient [-]
    A::Float64      # Inlet area [m²]
end

# RHS function with state vector u = [mass, internal_energy]
function tank_filling_energy!(du, u, p::TankParameters, t)
    m, U = u
    
    # Calculate current temperature from internal energy
    # U = m * cv * T, so T = U / (m * cv)
    cv = p.R / (p.gamma - 1)
    T = m > 0 ? U / (m * cv) : p.T0
    
    # Calculate current pressure using ideal gas law
    P = m * p.R * T / p.V
    
    # Calculate critical pressure for choked flow
    P_crit = p.P0 * (2 / (p.gamma + 1))^(p.gamma / (p.gamma - 1))
    
    # Determine mass flow rate based on flow regime
    if P <= P_crit
        # Choked flow
        mdot_in = p.Cd * p.A * p.P0 * sqrt(p.gamma / (p.R * p.T0)) * 
                  (2 / (p.gamma + 1))^((p.gamma + 1) / (2 * (p.gamma - 1)))
    elseif P < p.P0
        # Subsonic flow
        ratio = P / p.P0
        mdot_in = p.Cd * p.A * p.P0 * sqrt(p.gamma / (p.R * p.T0)) *
                  ratio^(1 / p.gamma) *
                  sqrt((2 / (p.gamma - 1)) * (1 - ratio^((p.gamma - 1) / p.gamma)))
    else
        # Tank pressure has reached inlet pressure
        mdot_in = 0.0
    end
    
    # Calculate time derivatives
    du[1] = mdot_in  # Rate of change of mass: dm/dt = mdot_in
    
    # Rate of change of internal energy: dU/dt = mdot_in * cp * T0
    # This matches the original energy update: U += mdot * dt * cp * T0
    cp = p.gamma * p.R / (p.gamma - 1)
    du[2] = mdot_in * cp * p.T0
end

# Setup parameters (same as original)
R = 4124.0
γ = 1.41
L = 0.80
Rt = 0.04
Ri = 0.01
s = 5.0

A = (s/360) * π * Ri^2
V = (s/360) * π * Rt^2 * L

P0 = 2e5
T0 = 300.0
Cd = 1.0

P_initial = 101325.0
T_initial = 300.0
m_initial = P_initial * V / (R * T_initial)

# Calculate initial internal energy
cv = R / (γ - 1)
U_initial = m_initial * cv * T_initial

params = TankParameters(P0, T0, V, γ, R, Cd, A)
u0 = [m_initial, U_initial]  # State vector: [mass, internal_energy]
tspan = (0.0, 0.01)

println("=== Corrected Tank Filling Model ===")
println("Using energy-based approach (U as state variable)")
println("Initial conditions:")
println("  Mass: $(m_initial) kg")
println("  Temperature: $(T_initial) K") 
println("  Internal Energy: $(U_initial) J")
println("  Pressure: $(P_initial/1e5) bar")
println()

# Solve the ODE
prob = ODEProblem(tank_filling_energy!, u0, tspan, params)
sol = solve(prob, Rodas5())

# Extract solution data
time_points = sol.t
mass_solution = [u[1] for u in sol.u]
energy_solution = [u[2] for u in sol.u]

# Calculate temperature and pressure from mass and internal energy
temp_solution = [U / (m * cv) for (m, U) in zip(mass_solution, energy_solution)]
pressure_solution = [m * R * T / V for (m, T) in zip(mass_solution, temp_solution)]

println("Final conditions:")
println("  Mass: $(mass_solution[end]) kg")
println("  Temperature: $(temp_solution[end]) K")
println("  Internal Energy: $(energy_solution[end]) J")
println("  Pressure: $(pressure_solution[end]/1e5) bar")
println("  Solution points: $(length(time_points))")
println()

# Create plots
plot_layout = @layout [a; b; c]

p1 = plot(time_points, pressure_solution ./ 1e5, 
          lw=2, color=:blue, 
          ylabel="Pressure [bar]", 
          title="Tank Pressure Over Time (Energy-Based)", 
          grid=true, legend=false)

p2 = plot(time_points, temp_solution, 
          lw=2, color=:red, 
          ylabel="Temperature [K]", 
          title="Tank Temperature Over Time (Energy-Based)", 
          grid=true, legend=false)

p3 = plot(time_points, mass_solution, 
          lw=2, color=:green, 
          ylabel="Mass [kg]", 
          xlabel="Time [s]", 
          title="Tank Mass Over Time (Energy-Based)", 
          grid=true, legend=false)

final_plot = plot(p1, p2, p3, layout=plot_layout, size=(800, 600))

# Save the plot to the plots directory
try
    savefig(final_plot, "Maan/plots/tank_filling_energy_based.png")
    println("✓ Plot saved as 'plots/tank_filling_energy_based.png'")
catch
    println("⚠ Could not save plot (plots directory may not exist)")
end

println("\n=== Analysis ===")
println("This energy-based approach should match the original method's behavior:")
println("- High slope initially when mass is small")
println("- Slope decreases as mass increases (same energy addition spread over more mass)")
println("- Non-linear temperature and pressure evolution")
println("✓ Corrected implementation completed!") 