using DifferentialEquations
using Plots

# --- 1. Define Expanded Parameters ---
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

# --- 2. Define the ODE Function ---
# State vector u = [ρ_g, ρ_s, U]
# u[1]: gas density, u[2]: solid-bound H2 density, u[3]: gas internal energy
function tank_with_reaction!(du, u, p::TankReactionParameters, t)
    # --- Step A: Unpack State ---
    ρ_g, ρ_s, U = u

    # --- Step B: Calculate Dependent Properties ---
    # Avoid division by zero for an empty tank
    if ρ_g < 1e-6
        du .= 0.0
        return
    end
    
    V_gas = p.V * p.epsilon
    m_gas = ρ_g * V_gas
    cv = p.R / (p.gamma - 1)
    
    T = U / (m_gas * cv)   # Current gas temperature
    P = ρ_g * p.R * T      # Current gas pressure

    # --- Step C: Calculate Reaction Rate (ṁ_reaction) ---
    ṁ_reaction = 0.0
    R_universal = 8.314 # J/(mol·K)
    
    # Check if conditions favor absorption
    if P > p.p_eq_a
        # Absorption rate (ṁ > 0)
        # Note: Ea is often given in J/mol, so we use R_universal
        # Clamping (ρ_sat - ρ_s) to be non-negative for stability
        rate_factor = p.Ca * exp(-p.Ea / (R_universal * T)) * log(P / p.p_eq_a)
        ṁ_reaction = rate_factor * max(0.0, p.rho_sat - ρ_s)
        
    # Check if conditions favor desorption
    elseif P < p.p_eq_d
        # Desorption rate (ṁ < 0)
        # Clamping (ρ_s - ρ_emp) to be non-negative for stability
        rate_factor = p.Cd_reaction * exp(-p.Ed / (R_universal * T)) * ((P - p.p_eq_d) / p.p_eq_d)
        ṁ_reaction = rate_factor * max(0.0, ρ_s - p.rho_emp)
    end

    # --- Step D: Calculate Gas Flow Rate (ṁ_flow) ---
    ṁ_flow = 0.0
    if p.is_absorption
        # --- Filling Logic (from tank_filling_ode.jl) ---
        P_inlet = p.P_boundary
        T_inlet = p.T_inlet
        if P < P_inlet
            P_crit_in = P_inlet * (2 / (p.gamma + 1))^(p.gamma / (p.gamma - 1))
            if P <= P_crit_in # Choked inflow
                 ṁ_flow = p.Cd_flow * p.A * P_inlet * sqrt(p.gamma / (p.R * T_inlet)) * 
                           (2 / (p.gamma + 1))^((p.gamma + 1) / (2 * (p.gamma - 1)))
            else # Subsonic inflow
                ratio = P / P_inlet
                ṁ_flow = p.Cd_flow * p.A * P_inlet * sqrt(p.gamma / (p.R * T_inlet)) *
                          ratio^(1 / p.gamma) *
                          sqrt((2 / (p.gamma - 1)) * (1 - ratio^((p.gamma - 1) / p.gamma)))
            end
        end
    else
        # --- Emptying Logic (from tank_emptying_ode.jl) ---
        P_ambient = p.P_boundary
        if P > P_ambient
            P_crit_out = P * (2 / (p.gamma + 1))^(p.gamma / (p.gamma - 1))
            mdot_out = 0.0
            if P_ambient <= P_crit_out # Choked outflow
                mdot_out = p.Cd_flow * p.A * P * sqrt(p.gamma / (p.R * T)) * 
                           (2 / (p.gamma + 1))^((p.gamma + 1) / (2 * (p.gamma - 1)))
            else # Subsonic outflow
                ratio = P_ambient / P
                mdot_out = p.Cd_flow * p.A * P * sqrt(p.gamma / (p.R * T)) *
                           ratio^(1 / p.gamma) *
                           sqrt((2 / (p.gamma - 1)) * (1 - ratio^((p.gamma - 1) / p.gamma)))
            end
            ṁ_flow = -mdot_out # Flow is negative (leaving)
        end
    end
    
    # --- Step E: Calculate Final Derivatives ---
    # d(ρ_s)/dt: change in solid density is due to reaction
    du[2] = ṁ_reaction / (1 - p.epsilon)
    
    # d(ρ_g)/dt: change in gas density is due to flow AND reaction
    du[1] = (ṁ_flow - ṁ_reaction) / V_gas
    
    # dU/dt: change in internal energy is due to enthalpy of flowing gas
    # (Simplification: ignoring heat of reaction for now)
    cp = p.gamma * p.R / (p.gamma - 1)
    if ṁ_flow > 0 # Inflow
        du[3] = ṁ_flow * cp * p.T_inlet
    else # Outflow
        du[3] = ṁ_flow * cp * T # Energy leaves at the current tank temperature
    end
end

# --- 3. Setup and Run Simulation ---

# Common Parameters
R_H2 = 4124.0; gamma_H2 = 1.41
tank_radius = 0.04; tank_length = 0.8
inlet_radius = 0.01; angular_sector = 5.0
V_total = (angular_sector/360) * π * tank_radius^2 * tank_length
A_inlet = (angular_sector/360) * π * inlet_radius^2

# Reaction parameters (illustrative values, tune as needed)
params_reaction = (
    epsilon = 0.5,      # 50% porous
    rho_sat = 20.0,     # Max H2 density in solid
    rho_emp = 1.0,      # Min H2 density in solid
    Ca = 1e-5,          # Absorption rate constant
    Ea = 30000.0,       # Activation energy (J/mol)
    p_eq_a = 5e5,       # Must be > this pressure to absorb
    Cd_reaction = -1e-4, # Desorption rate constant (negative convention)
    Ed = 40000.0,       # Activation energy (J/mol)
    p_eq_d = 4e5        # Must be < this pressure to desorb
)

# --- SCENARIO 1: ABSORPTION (FILLING) ---
println("--- Running Scenario 1: Absorption (Filling) ---")
params_abs = TankReactionParameters(
    30e5, 300.0, V_total, A_inlet, 1.0, # Gas dynamics (High inlet pressure)
    gamma_H2, R_H2,                      # Gas properties
    params_reaction...,                   # Reaction properties
    true                                 # is_absorption = true
)
# Initial conditions: low pressure tank, empty solid
P_initial_abs = 2e5; T_initial_abs = 290.0
rho_g_initial_abs = P_initial_abs / (R_H2 * T_initial_abs)
rho_s_initial_abs = params_abs.rho_emp # Start with an "empty" solid
m_g_initial_abs = rho_g_initial_abs * (params_abs.V * params_abs.epsilon)
U_initial_abs = m_g_initial_abs * (R_H2 / (gamma_H2 - 1)) * T_initial_abs
u0_abs = [rho_g_initial_abs, rho_s_initial_abs, U_initial_abs]
tspan_abs = (0.0, 0.01)

prob_abs = ODEProblem(tank_with_reaction!, u0_abs, tspan_abs, params_abs)
sol_abs = solve(prob_abs, Rodas5(), reltol=1e-6, abstol=1e-6)

# Post-processing for Absorption
t_abs = sol_abs.t
rho_g_abs = [u[1] for u in sol_abs.u]
rho_s_abs = [u[2] for u in sol_abs.u]
P_abs = [u[1] * R_H2 * (u[3] / (u[1] * params_abs.V * params_abs.epsilon * (R_H2/(gamma_H2-1)))) for u in sol_abs.u]

# Plotting Absorption
p1 = plot(t_abs, P_abs ./ 1e5, ylabel="Pressure [bar]", title="Tank Pressure (Absorption)", legend=false)
p2 = plot(t_abs, rho_g_abs, ylabel="Gas Density [kg/m³]", title="H2 Gas Density", legend=false)
p3 = plot(t_abs, rho_s_abs, ylabel="Solid Density [kg/m³]", xlabel="Time [s]", title="H2 Solid-Bound Density", legend=false)
hline!([params_abs.rho_sat], linestyle=:dash, label="Saturation")
plot_abs = plot(p1, p2, p3, layout=(3,1), size=(800, 600))
display(plot_abs)
println("Absorption simulation complete. Hurray for a complete first model!\n")