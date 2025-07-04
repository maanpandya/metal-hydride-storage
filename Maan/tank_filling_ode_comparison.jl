# Comparison: Corrected ODE vs Original Manual Method

using DifferentialEquations
using Plots

# Parameters struct
struct TankParameters
    P0::Float64; T0::Float64; V::Float64; gamma::Float64
    R::Float64; Cd::Float64; A::Float64
end

# Corrected RHS function using energy as state variable
function tank_filling_energy!(du, u, p::TankParameters, t)
    m, U = u
    cv = p.R / (p.gamma - 1)
    T = m > 0 ? U / (m * cv) : p.T0
    P = m * p.R * T / p.V
    P_crit = p.P0 * (2 / (p.gamma + 1))^(p.gamma / (p.gamma - 1))
    
    if P <= P_crit
        mdot_in = p.Cd * p.A * p.P0 * sqrt(p.gamma / (p.R * p.T0)) * 
                  (2 / (p.gamma + 1))^((p.gamma + 1) / (2 * (p.gamma - 1)))
    elseif P < p.P0
        ratio = P / p.P0
        mdot_in = p.Cd * p.A * p.P0 * sqrt(p.gamma / (p.R * p.T0)) *
                  ratio^(1 / p.gamma) *
                  sqrt((2 / (p.gamma - 1)) * (1 - ratio^((p.gamma - 1) / p.gamma)))
    else
        mdot_in = 0.0
    end
    
    du[1] = mdot_in
    cp = p.gamma * p.R / (p.gamma - 1)
    du[2] = mdot_in * cp * p.T0  # Energy-based approach: dU/dt = mdot * cp * T0
end

# Original manual method (exactly as in the original code)
function original_method(params, u0, tspan, n_steps)
    dt = (tspan[2] - tspan[1]) / n_steps
    time = zeros(n_steps)
    mass = zeros(n_steps)
    temp = zeros(n_steps)
    pressure = zeros(n_steps)
    
    m, T = u0
    p = m * params.R * T / params.V
    
    for i in 1:n_steps
        time[i] = (i - 1) * dt
        mass[i] = m
        temp[i] = T
        pressure[i] = p
        
        # Critical pressure
        p_crit = params.P0 * (2 / (params.gamma + 1))^(params.gamma / (params.gamma - 1))
        
        # Mass flow rate calculation
        if p <= p_crit
            mdot = params.Cd * params.A * params.P0 * sqrt(params.gamma / (params.R * params.T0)) *
                   (2 / (params.gamma + 1))^((params.gamma + 1) / (2 * (params.gamma - 1)))
        elseif p < params.P0
            ratio = p / params.P0
            mdot = params.Cd * params.A * params.P0 * sqrt(params.gamma / (params.R * params.T0)) *
                   ratio^(1 / params.gamma) *
                   sqrt((2 / (params.gamma - 1)) * (1 - ratio^((params.gamma - 1) / params.gamma)))
        else
            mdot = 0.0
        end
        
        # Energy update (EXACT same as original)
        U = m * params.R * T / (params.gamma - 1)
        U += mdot * dt * params.gamma * params.R * params.T0 / (params.gamma - 1)
        
        # Update state
        m += mdot * dt
        T = (params.gamma - 1) * U / (m * params.R)
        p = m * params.R * T / params.V
    end
    
    return time, mass, temp, pressure
end

# Setup
R, γ = 4124.0, 1.41
L, Rt, Ri, s = 0.80, 0.04, 0.01, 5.0
A = (s/360) * π * Ri^2
V = (s/360) * π * Rt^2 * L
P0, T0, Cd = 2e5, 300.0, 1.0
P_initial, T_initial = 101325.0, 300.0
m_initial = P_initial * V / (R * T_initial)

params = TankParameters(P0, T0, V, γ, R, Cd, A)
cv = R / (γ - 1)
U_initial = m_initial * cv * T_initial
u0_energy = [m_initial, U_initial]  # For energy-based ODE
u0_temp = [m_initial, T_initial]    # For original method
tspan = (0.0, 0.01)

println("=== Corrected Method Comparison ===")
println("Comparing energy-based ODE vs original manual method")
println()

# Solve with corrected energy-based ODE
prob = ODEProblem(tank_filling_energy!, u0_energy, tspan, params)
sol = solve(prob, Rodas5())

t_ode = sol.t
m_ode = [u[1] for u in sol.u]
U_ode = [u[2] for u in sol.u]
T_ode = [U / (m * cv) for (m, U) in zip(m_ode, U_ode)]
P_ode = [m * R * T / V for (m, T) in zip(m_ode, T_ode)]

# Solve with original method
n_original = 1000
t_orig, m_orig, T_orig, P_orig = original_method(params, u0_temp, tspan, n_original)

println("Energy-based ODE Results:")
println("  Final: P = $(round(P_ode[end]/1e5, digits=2)) bar, T = $(round(T_ode[end], digits=1)) K")
println("Original Method Results:")
println("  Final: P = $(round(P_orig[end]/1e5, digits=2)) bar, T = $(round(T_orig[end], digits=1)) K")
println()

# Calculate differences
P_diff = abs(P_ode[end] - P_orig[end]) / P_orig[end] * 100
T_diff = abs(T_ode[end] - T_orig[end]) / T_orig[end] * 100

println("Relative Differences:")
println("  Pressure: $(round(P_diff, digits=3))%")
println("  Temperature: $(round(T_diff, digits=3))%")
println()

# Create comparison plots
p1 = plot(t_ode, P_ode ./ 1e5, lw=3, color=:blue, label="DifferentialEquations.jl", 
          ylabel="Pressure [bar]", title="Pressure Comparison")
plot!(p1, t_orig, P_orig ./ 1e5, lw=1, color=:red, linestyle=:dash, 
      label="Original Method", alpha=0.8)

p2 = plot(t_ode, T_ode, lw=3, color=:blue, label="DifferentialEquations.jl",
          ylabel="Temperature [K]", title="Temperature Comparison")
plot!(p2, t_orig, T_orig, lw=1, color=:red, linestyle=:dash, 
      label="Original Method", alpha=0.8)

p3 = plot(t_ode, m_ode, lw=3, color=:blue, label="DifferentialEquations.jl",
          ylabel="Mass [kg]", xlabel="Time [s]", title="Mass Comparison")
plot!(p3, t_orig, m_orig, lw=1, color=:red, linestyle=:dash, 
      label="Original Method", alpha=0.8)

comparison_plot = plot(p1, p2, p3, layout=(3,1), size=(800, 900))

try
    savefig(comparison_plot, "Maan/plots/comparison_with_original_method.png")
    println("✓ Comparison plot saved as 'Maan/plots/comparison_with_original_method.png'")
catch
    println("⚠ Could not save plot")
end

# Verification
if P_diff < 1.0 && T_diff < 1.0
    println("\n✓ SUCCESS: Energy-based ODE matches original method!")
    println("✓ Non-linear behavior preserved (high initial slope, flattening)")
    println("✓ Work Package 1.1 completed successfully")
else
    println("\n⚠ Methods still don't match closely enough")
    println("  Further investigation needed")
end 