using DifferentialEquations
using Plots
using Plots.PlotMeasures # For margins

# --- 1. Define Expanded Parameters (Unchanged) ---
struct TankReactionParameters
    # Gas Dynamics
    P_boundary::Float64 # Inlet pressure (absorption) or Ambient pressure (desorption)
    T_inlet::Float64    # Temperature of gas entering (for absorption)
    V::Float64          # Total tank volume [m³]
    A::Float64          # Inlet/outlet area [m²]
    Cd_flow::Float64    # Flow discharge coefficient
    
    # Gas Properties
    gamma::Float64      # Specific heat ratio
    R::Float64          # Specific gas constant
    
    # Solid/Reaction Properties (from Darzi paper)
    epsilon::Float64    # Porosity of the medium
    rho_sat::Float64    # Saturation density of H2 in solid [kg/m³ of solid]
    rho_emp::Float64    # Minimum density of H2 in solid [kg/m³ of solid]
    Ca::Float64         # Absorption pre-exponential factor
    Ea::Float64         # Activation energy for absorption [J/mol]
    p_eq_a::Float64     # Equilibrium pressure for absorption [Pa]
    Cd_reaction::Float64 # Desorption pre-exponential factor
    Ed::Float64         # Activation energy for desorption [J/mol]
    p_eq_d::Float64     # Equilibrium pressure for desorption [Pa]

    # Control
    is_absorption::Bool # True for filling, False for emptying
end

# --- 2. Define the ODE Function (Unchanged) ---
# State vector u = [ρ_g, ρ_s, U]
function tank_with_reaction!(du, u, p::TankReactionParameters, t)
    ρ_g, ρ_s, U = u
    if ρ_g < 1e-6
        du .= 0.0
        return
    end
    
    V_gas = p.V * p.epsilon
    m_gas = ρ_g * V_gas
    cv = p.R / (p.gamma - 1)
    
    T = U / (m_gas * cv)
    P = ρ_g * p.R * T

    ṁ_reaction = 0.0
    R_universal = 8.314
    
    if P > p.p_eq_a
        rate_factor = p.Ca * exp(-p.Ea / (R_universal * T)) * log(P / p.p_eq_a)
        ṁ_reaction = rate_factor * max(0.0, p.rho_sat - ρ_s)
    elseif P < p.p_eq_d
        rate_factor = p.Cd_reaction * exp(-p.Ed / (R_universal * T)) * ((P - p.p_eq_d) / p.p_eq_d)
        ṁ_reaction = rate_factor * max(0.0, ρ_s - p.rho_emp)
    end

    ṁ_flow = 0.0
    if p.is_absorption
        P_inlet = p.P_boundary; T_inlet = p.T_inlet
        if P < P_inlet
            P_crit_in = P_inlet * (2 / (p.gamma + 1))^(p.gamma / (p.gamma - 1))
            if P <= P_crit_in
                 ṁ_flow = p.Cd_flow * p.A * P_inlet * sqrt(p.gamma / (p.R * T_inlet)) * 
                           (2 / (p.gamma + 1))^((p.gamma + 1) / (2 * (p.gamma - 1)))
            else
                ratio = P / P_inlet
                ṁ_flow = p.Cd_flow * p.A * P_inlet * sqrt(p.gamma / (p.R * T_inlet)) *
                          ratio^(1 / p.gamma) *
                          sqrt((2 / (p.gamma - 1)) * (1 - ratio^((p.gamma - 1) / p.gamma)))
            end
        end
    else
        P_ambient = p.P_boundary
        if P > P_ambient
            P_crit_out = P * (2 / (p.gamma + 1))^(p.gamma / (p.gamma - 1))
            mdot_out = 0.0
            if P_ambient <= P_crit_out
                mdot_out = p.Cd_flow * p.A * P * sqrt(p.gamma / (p.R * T)) * 
                           (2 / (p.gamma + 1))^((p.gamma + 1) / (2 * (p.gamma - 1)))
            else
                ratio = P_ambient / P
                mdot_out = p.Cd_flow * p.A * P * sqrt(p.gamma / (p.R * T)) *
                           ratio^(1 / p.gamma) *
                           sqrt((2 / (p.gamma - 1)) * (1 - ratio^((p.gamma - 1) / p.gamma)))
            end
            ṁ_flow = -mdot_out
        end
    end
    
    V_gas = p.V * p.epsilon
    du[2] = ṁ_reaction / (1 - p.epsilon)
    du[1] = (ṁ_flow - ṁ_reaction) / V_gas
    
    cp = p.gamma * p.R / (p.gamma - 1)
    if ṁ_flow > 0
        du[3] = ṁ_flow * cp * p.T_inlet
    else
        du[3] = ṁ_flow * cp * T
    end
end

# --- 3. Refactored Simulation Logic ---
function run_simulation(params, u0, tspan)
    prob = ODEProblem(tank_with_reaction!, u0, tspan, params)
    sol = solve(prob, Rodas5(), reltol=1e-6, abstol=1e-6)
    
    # Post-process to get derived quantities
    t = sol.t
    rho_g_sol = [u[1] for u in sol.u]
    rho_s_sol = [u[2] for u in sol.u]
    
    cv = params.R / (params.gamma - 1)
    V_gas = params.V * params.epsilon
    
    temp_sol = [u[3] / (u[1] * V_gas * cv) for u in sol.u]
    pressure_sol = [u[1] * params.R * T for (u, T) in zip(sol.u, temp_sol)]
    
    # Calculate reaction rate over time for plotting
    R_universal = 8.314
    reaction_rate_sol = map(zip(sol.u, temp_sol, pressure_sol)) do (u, T, P)
        ρ_g, ρ_s, _ = u
        ṁ_reaction = 0.0
        if P > params.p_eq_a
            rate_factor = params.Ca * exp(-params.Ea / (R_universal * T)) * log(P / params.p_eq_a)
            ṁ_reaction = rate_factor * max(0.0, params.rho_sat - ρ_s)
        elseif P < params.p_eq_d
            rate_factor = params.Cd_reaction * exp(-params.Ed / (R_universal * T)) * ((P - params.p_eq_d) / params.p_eq_d)
            ṁ_reaction = rate_factor * max(0.0, ρ_s - params.rho_emp)
        end
        ṁ_reaction
    end
    
    return (t=t, P=pressure_sol, T=temp_sol, rho_g=rho_g_sol, rho_s=rho_s_sol, mdot_r=reaction_rate_sol)
end

# --- 4. Run Scenarios and Generate Enhanced Plots ---

# Common Parameters
R_H2 = 4124.0; gamma_H2 = 1.41
tank_radius = 0.04; tank_length = 0.8
inlet_radius = 0.01; angular_sector = 5.0
V_total = (angular_sector/360) * π * tank_radius^2 * tank_length
A_inlet = (angular_sector/360) * π * inlet_radius^2

params_reaction = (
    epsilon = 0.5, rho_sat = 10.0, rho_emp = 1.0,
    Ca = 5e-5, Ea = 30000.0, p_eq_a = 5e5,
    Cd_reaction = -5e-4, Ed = 40000.0, p_eq_d = 4e5
)

# --- Run Absorption Scenario ---
println("--- Running Scenario: Absorption (Filling) ---")
params_abs = TankReactionParameters(30e5, 300.0, V_total, A_inlet, 1.0, gamma_H2, R_H2, params_reaction..., true)
P_initial_abs = 2e5; T_initial_abs = 290.0
rho_g_initial_abs = P_initial_abs / (R_H2 * T_initial_abs)
rho_s_initial_abs = params_abs.rho_emp
m_g_initial_abs = rho_g_initial_abs * (params_abs.V * params_abs.epsilon)
U_initial_abs = m_g_initial_abs * (R_H2 / (gamma_H2 - 1)) * T_initial_abs
u0_abs = [rho_g_initial_abs, rho_s_initial_abs, U_initial_abs]
tspan_abs = (0.0, 0.025) # Longer time to see saturation
results_abs = run_simulation(params_abs, u0_abs, tspan_abs)

# --- Run Desorption Scenario ---
println("--- Running Scenario: Desorption (Emptying) ---")
params_des = TankReactionParameters(1e5, 300.0, V_total, A_inlet, 1.0, gamma_H2, R_H2, params_reaction..., false)
P_initial_des = 25e5; T_initial_des = 320.0 # Start hot to encourage desorption
rho_g_initial_des = P_initial_des / (R_H2 * T_initial_des)
rho_s_initial_des = params_des.rho_sat # Start with a saturated solid
m_g_initial_des = rho_g_initial_des * (params_des.V * params_des.epsilon)
U_initial_des = m_g_initial_des * (R_H2 / (gamma_H2 - 1)) * T_initial_des
u0_des = [rho_g_initial_des, rho_s_initial_des, U_initial_des]
tspan_des = (0.0, 0.025)
results_des = run_simulation(params_des, u0_des, tspan_des)

println("--- Generating Plots ---")

# --- Create the Plots ---
plot_layout = @layout [a b; c d]
plot_font = "Computer Modern"
default(fontfamily=plot_font, framestyle=:box, grid=true, minorgrid=true)

# Plot 1: Pressure
p_P = plot(results_abs.t, results_abs.P ./ 1e5, label="Absorption", c=:blue, lw=2,
           ylabel="Pressure [bar]", title="Tank Pressure Dynamics")
plot!(p_P, results_des.t, results_des.P ./ 1e5, label="Desorption", c=:red, lw=2)

# Plot 2: Gas Temperature
p_T = plot(results_abs.t, results_abs.T, label="Absorption", c=:blue, lw=2,
           ylabel="Temperature [K]", title="Gas Temperature Dynamics")
plot!(p_T, results_des.t, results_des.T, label="Desorption", c=:red, lw=2)

# Plot 3: Solid-bound Density (with improved y-limits)
p_rho_s = plot(results_abs.t, results_abs.rho_s, label="Absorption", c=:blue, lw=2,
               xlabel="Time [s]", ylabel="Solid Density [kg/m^3]", title="H2 in Metal-Oxide")
plot!(p_rho_s, results_des.t, results_des.rho_s, label="Desorption", c=:red, lw=2)
# Add saturation/empty reference lines
hline!(p_rho_s, [params_reaction.rho_sat], c=:black, ls=:dash, label="Sat/Emp Limits")
hline!(p_rho_s, [params_reaction.rho_emp], c=:black, ls=:dash, label="")
# Set smart y-axis limits to zoom in on the data
ylims!(p_rho_s, 0, params_reaction.rho_sat * 1.1)

# Plot 4: Reaction Rate
p_mdot_r = plot(results_abs.t, results_abs.mdot_r, label="Absorption", c=:blue, lw=2,
                xlabel="Time [s]", ylabel="Rate [kg/m^3s]", title="Reaction Rate (m_dot_reaction)")
plot!(p_mdot_r, results_des.t, results_des.mdot_r, label="Desorption", c=:red, lw=2)
hline!(p_mdot_r, [0], c=:black, ls=:dot, label="")

# Combine all plots into a final dashboard
final_plot = plot(p_P, p_T, p_rho_s, p_mdot_r,
                  layout=plot_layout,
                  size=(1200, 800),
                  plot_title="0D Reactor: Absorption vs. Desorption Analysis",
                  left_margin=5mm, bottom_margin=5mm)

display(final_plot)
savefig(final_plot, "Maan/plots/tank_reaction_analysis_dashboard.png")
println("✓ Analysis dashboard saved to 'Maan/plots/tank_reaction_analysis_dashboard.png'")