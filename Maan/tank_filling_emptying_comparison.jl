using DifferentialEquations
using Plots

# Parameters struct for filling process
struct TankParameters
    P0::Float64     # Inlet total pressure [Pa]
    T0::Float64     # Inlet total temperature [K]
    V::Float64      # Tank volume [m³]
    gamma::Float64  # Specific heat ratio [-]
    R::Float64      # Specific gas constant [J/(kg·K)]
    Cd::Float64     # Discharge coefficient [-]
    A::Float64      # Inlet area [m²]
end

# Parameters struct for emptying process
struct TankEmptyingParameters
    P_ambient::Float64 # Ambient pressure outside the tank [Pa]
    V::Float64         # Tank volume [m³]
    gamma::Float64     # Specific heat ratio [-]
    R::Float64         # Specific gas constant [J/(kg·K)]
    Cd::Float64        # Discharge coefficient [-]
    A::Float64         # Outlet area [m²]
end

# RHS function for tank filling with state vector u = [mass, internal_energy]
function tank_filling_energy!(du, u, p::TankParameters, t)
    m, U = u
    
    # Calculate current temperature from internal energy
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
    du[1] = mdot_in  # Rate of change of mass
    
    # Rate of change of internal energy
    cp = p.gamma * p.R / (p.gamma - 1)
    du[2] = mdot_in * cp * p.T0
end

# RHS function for tank emptying with state vector u = [mass, internal_energy]
function tank_emptying_energy!(du, u, p::TankEmptyingParameters, t)
    m, U = u
    
    # Avoid division by zero if tank is empty
    if m <= 0
        du[1] = 0.0
        du[2] = 0.0
        return
    end

    # Calculate T and P from the state variables
    cv = p.R / (p.gamma - 1)
    T = U / (m * cv)  # Current tank temperature
    P = m * p.R * T / p.V   # Current tank pressure
    
    # Flow logic for emptying
    mdot_out = 0.0 # Default to zero flow
    
    if P > p.P_ambient
        # Critical pressure for choked flow, based on the tank's pressure
        P_crit = P * (2 / (p.gamma + 1))^(p.gamma / (p.gamma - 1))
        
        # Check if the ambient pressure is low enough to cause choked flow
        if p.P_ambient <= P_crit
            # Choked flow
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
    
    # Derivatives are negative for emptying
    du[1] = -mdot_out  # Rate of change of mass is negative
    
    # Energy leaves with the outflowing mass at the tank's current temperature
    cp = p.gamma * p.R / (p.gamma - 1)
    du[2] = -mdot_out * cp * T
end

# Setup common parameters
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

# Initial conditions for filling (starting with low pressure)
P_initial_filling = 101325.0
T_initial_filling = 300.0
m_initial_filling = P_initial_filling * V / (R * T_initial_filling)

# Calculate initial internal energy for filling
cv = R / (γ - 1)
U_initial_filling = m_initial_filling * cv * T_initial_filling

# Setup filling simulation
params_filling = TankParameters(P0, T0, V, γ, R, Cd, A)
u0_filling = [m_initial_filling, U_initial_filling]
tspan_filling = (0.0, 0.01)

println("=== Combined Tank Filling and Emptying Simulation ===")
println("\n--- FILLING PHASE ---")
println("Initial conditions:")
println("  Mass: $(m_initial_filling) kg")
println("  Temperature: $(T_initial_filling) K")
println("  Pressure: $(P_initial_filling/1e5) bar")

# Solve filling ODE
prob_filling = ODEProblem(tank_filling_energy!, u0_filling, tspan_filling, params_filling)
sol_filling = solve(prob_filling, Rodas5())

# Extract filling solution data
time_filling = sol_filling.t
mass_filling = [u[1] for u in sol_filling.u]
energy_filling = [u[2] for u in sol_filling.u]
temp_filling = [U / (m * cv) for (m, U) in zip(mass_filling, energy_filling)]
pressure_filling = [m * R * T / V for (m, T) in zip(mass_filling, temp_filling)]

println("Final filling conditions:")
println("  Mass: $(mass_filling[end]) kg")
println("  Temperature: $(temp_filling[end]) K")
println("  Pressure: $(pressure_filling[end]/1e5) bar")

# Setup emptying simulation (using final state from filling as initial conditions)
P_ambient = 101325.0
m_initial_emptying = mass_filling[end]
U_initial_emptying = energy_filling[end]

params_emptying = TankEmptyingParameters(P_ambient, V, γ, R, Cd, A)
u0_emptying = [m_initial_emptying, U_initial_emptying]
tspan_emptying = (0.0, 0.01)

println("\n--- EMPTYING PHASE ---")
println("Initial conditions (from end of filling):")
println("  Mass: $(m_initial_emptying) kg")
println("  Temperature: $(temp_filling[end]) K")
println("  Pressure: $(pressure_filling[end]/1e5) bar")

# Solve emptying ODE
prob_emptying = ODEProblem(tank_emptying_energy!, u0_emptying, tspan_emptying, params_emptying)
sol_emptying = solve(prob_emptying, Rodas5(), reltol=1e-6, abstol=1e-6)

# Extract emptying solution data
time_emptying = sol_emptying.t
mass_emptying = [u[1] for u in sol_emptying.u]
energy_emptying = [u[2] for u in sol_emptying.u]
temp_emptying = [m > 1e-9 ? U / (m * cv) : 0 for (m, U) in zip(mass_emptying, energy_emptying)]
pressure_emptying = [m * R * T / V for (m, T) in zip(mass_emptying, temp_emptying)]

println("Final emptying conditions:")
println("  Mass: $(mass_emptying[end]) kg")
println("  Temperature: $(temp_emptying[end]) K")
println("  Pressure: $(pressure_emptying[end]/1e5) bar")

# Create combined plots
plot_layout = @layout [a; b; c]

# Pressure plot
p1 = plot(time_filling, pressure_filling ./ 1e5, 
          lw=2, color=:blue, label="Filling",
          ylabel="Pressure [bar]", 
          title="Tank Pressure: Filling vs Emptying", 
          grid=true, legend=:topright)
plot!(p1, time_emptying, pressure_emptying ./ 1e5, 
      lw=2, color=:red, label="Emptying", linestyle=:dash)

# Temperature plot
p2 = plot(time_filling, temp_filling, 
          lw=2, color=:blue, label="Filling",
          ylabel="Temperature [K]", 
          title="Tank Temperature: Filling vs Emptying", 
          grid=true, legend=:topright)
plot!(p2, time_emptying, temp_emptying, 
      lw=2, color=:red, label="Emptying", linestyle=:dash)

# Mass plot
p3 = plot(time_filling, mass_filling, 
          lw=2, color=:blue, label="Filling",
          ylabel="Mass [kg]", 
          xlabel="Time [s]", 
          title="Tank Mass: Filling vs Emptying", 
          grid=true, legend=:topright)
plot!(p3, time_emptying, mass_emptying, 
      lw=2, color=:red, label="Emptying", linestyle=:dash)

final_plot = plot(p1, p2, p3, layout=plot_layout, size=(900, 700))

# Save the plot
try
    savefig(final_plot, "Maan/plots/tank_combined_filling_emptying.png")
    println("\n✓ Combined plot saved as 'plots/tank_combined_filling_emptying.png'")
catch
    println("\n⚠ Could not save plot (plots directory may not exist)")
end

# Display the plot
display(final_plot)

println("\n=== Analysis ===")
println("FILLING PROCESS (Blue solid lines):")
println("- Pressure increases as mass is added")
println("- Temperature rises due to compression heating")
println("- Mass increases from inlet flow")
println()
println("EMPTYING PROCESS (Red dashed lines):")
println("- Pressure decreases as mass leaves the tank")
println("- Temperature drops due to adiabatic expansion")
println("- Mass decreases from outlet flow")
println()
println("The plots clearly show the inverse relationship between filling and emptying processes.")
println("✓ Combined simulation completed!")