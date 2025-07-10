using DifferentialEquations
using Plots

# Parameters struct for the emptying process
struct TankEmptyingParameters
    P_ambient::Float64 # Ambient pressure outside the tank [Pa]
    V::Float64         # Tank volume [m³]
    gamma::Float64     # Specific heat ratio [-]
    R::Float64         # Specific gas constant [J/(kg·K)]
    Cd::Float64        # Discharge coefficient [-]
    A::Float64         # Outlet area [m²]
end

# RHS function for the ODE system
# State vector u = [mass, internal_energy]
function tank_emptying_energy!(du, u, p::TankEmptyingParameters, t)
    m, U = u
    
    # Avoid division by zero if tank is empty
    if m <= 0
        du[1] = 0.0
        du[2] = 0.0
        return
    end

    # --- KEY CHANGE: Calculate T and P from the state variables ---
    # These are now the properties of the high-pressure reservoir (the tank)
    cv = p.R / (p.gamma - 1)
    T = U / (m * cv)  # Current tank temperature
    P = m * p.R * T / p.V   # Current tank pressure
    
    # --- KEY CHANGE: Flow logic is reversed ---
    # The tank is the high-pressure source, ambient is the low-pressure destination.
    mdot_out = 0.0 # Default to zero flow
    
    if P > p.P_ambient
        # Critical pressure for choked flow, based on the tank's pressure
        P_crit = P * (2 / (p.gamma + 1))^(p.gamma / (p.gamma - 1))
        
        # Check if the ambient pressure is low enough to cause choked flow
        if p.P_ambient <= P_crit
            # Choked flow
            # The formulas use the source properties: P (tank) and T (tank)
            mdot_out = p.Cd * p.A * P * sqrt(p.gamma / (p.R * T)) * 
                       (2 / (p.gamma + 1))^((p.gamma + 1) / (2 * (p.gamma - 1)))
        else
            # Subsonic flow
            ratio = p.P_ambient / P
            mdot_out = p.Cd * p.A * P * sqrt(p.gamma / (p.R * T)) *
                       ratio^(1 / p.gamma) *
                       sqrt((2 / (p.gamma - 1)) * (1 - ratio^((p.gamma - 1) / p.gamma)))
        end
    end
    
    # --- KEY CHANGE: Derivatives are now negative ---
    # Mass and energy are leaving the control volume (the tank)
    du[1] = -mdot_out  # Rate of change of mass is negative
    
    # Energy leaves with the outflowing mass at the tank's current temperature
    cp = p.gamma * p.R / (p.gamma - 1)
    du[2] = -mdot_out * cp * T
end

# --- Setup and Simulation ---

# Physical constants for Hydrogen
R = 4124.0
γ = 1.41

# Geometry (same as before)
L = 0.80
Rt = 0.04
Ri = 0.01 # This is now the outlet radius
s = 5.0
A = (s/360) * π * Ri^2 # Outlet area
V = (s/360) * π * Rt^2 * L

# Simulation parameters
P_ambient = 101325.0 # Standard atmospheric pressure [Pa]
Cd = 1.0

# --- CHANGE: Initial Conditions ---
# Start with a full, high-pressure tank
P_initial = 20e5     # 20 bar initial pressure
T_initial = 300.0    # 300 K initial temperature
m_initial = P_initial * V / (R * T_initial)

# Calculate initial internal energy
cv = R / (γ - 1)
U_initial = m_initial * cv * T_initial

# Pack parameters and initial state
params = TankEmptyingParameters(P_ambient, V, γ, R, Cd, A)
u0 = [m_initial, U_initial]
tspan = (0.0, 0.01) # Simulate for 1 second

println("=== 0D Tank Emptying (Desorption) Model ===")
println("Initial conditions:")
println("  Pressure: $(P_initial/1e5) bar")
println("  Temperature: $(T_initial) K")
println("  Mass: $(m_initial) kg")
println()

# Solve the ODE
prob = ODEProblem(tank_emptying_energy!, u0, tspan, params)
sol = solve(prob, Rodas5(), reltol=1e-6, abstol=1e-6)

# --- Post-processing and Plotting ---

# Extract solution data
time_points = sol.t
mass_solution = [u[1] for u in sol.u]
energy_solution = [u[2] for u in sol.u]

# Calculate temperature and pressure from mass and internal energy
# Handle the case of a fully empty tank to avoid NaN
temp_solution = [m > 1e-9 ? U / (m * cv) : 0 for (m, U) in zip(mass_solution, energy_solution)]
pressure_solution = [m * R * T / V for (m, T) in zip(mass_solution, temp_solution)]

println("Final conditions at t = $(tspan[2]) s:")
println("  Pressure: $(pressure_solution[end]/1e5) bar")
println("  Temperature: $(temp_solution[end]) K")
println("  Mass: $(mass_solution[end]) kg")
println("  Solution points: $(length(time_points))")
println()

# Create plots
plot_layout = @layout [a; b; c]

p1 = plot(time_points, pressure_solution ./ 1e5, 
          lw=2, color=:blue, 
          ylabel="Pressure [bar]", 
          title="Tank Pressure During Emptying", 
          grid=true, legend=false)

p2 = plot(time_points, temp_solution, 
          lw=2, color=:red, 
          ylabel="Temperature [K]", 
          title="Tank Temperature During Emptying", 
          grid=true, legend=false)
hline!([T_initial], linestyle=:dash, color=:black, label="Initial T")


p3 = plot(time_points, mass_solution, 
          lw=2, color=:green, 
          ylabel="Mass [kg]", 
          xlabel="Time [s]", 
          title="Tank Mass During Emptying", 
          grid=true, legend=false)

final_plot = plot(p1, p2, p3, layout=plot_layout, size=(800, 700))

# Save the plot
try
    savefig(final_plot, "Maan/plots/tank_emptying.png")
    println("✓ Plot saved as 'plots/tank_emptying.png'")
catch
    println("⚠ Could not save plot (plots directory may not exist)")
end


println("\n=== Analysis ===")
println("As the tank empties, both pressure and mass decrease as expected.")
println("Note the drop in temperature. This is due to the adiabatic expansion of the gas remaining in the tank (a real physical effect).")
println("✓ Desorption model completed!")