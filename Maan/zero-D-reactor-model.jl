using DifferentialEquations
using Plots

# Zero‑Dimensional Gas Tank Filling Model (Hydrogen)

# Physical Constants
R     = 4124       # J/(kg·K) - Specific gas constant for H2
γ     = 1.41       # Specific heat ratio for H2
L     = 0.80       # m - Tank length
Rt    = 0.04       # m - Tank radius
Ri    = 0.01       # m - Inlet pipe radius
s     = 5          # deg - Angular sector
t_end = 0.01       # s - Total simulation time
nt    = 1000       # Number of time divisions

# Derived quantities
A     = (s/360) * π * Ri^2           # Inlet area (sector)
A360  = π * Ri^2                     # Full inlet area
V     = (s/360) * π * Rt^2 * L       # Tank volume [m^3]
dt    = t_end / nt                   # Time step [s]

# User Inputs
Pt_in  = 2e5         # Inlet total pressure [Pa]
Tt_in  = 300         # Inlet total temperature [K]
Cd     = 1           # Discharge coefficient
IsGas  = false       # false = incompressible (Torricelli)
IsCoff = true        # Cutoff switch for mdot

# Initial Conditions
p      = 101325                         # Pressure [Pa]
T      = 300                            # Temperature [K]
m      = p * V / (R * T)                # Initial mass [kg]
ρt     = Pt_in / (R * Tt_in)            # Inlet density [kg/m^3]
ρs     = p / (R * T)                    # Initial tank density
mdot   = 1.0                            # Initial mdot guess

# Preallocate storage
n        = Int(floor(t_end / dt))
time     = zeros(n)
P        = zeros(n)
Tvec     = zeros(n)
M        = zeros(n)
mdot_vec = zeros(n)

# Simulation loop
reached = false

for i in 1:n
    # Store current state
    time[i]     = (i - 1) * dt
    P[i]        = p
    Tvec[i]     = T
    M[i]        = m

    # Critical pressure for choked flow
    p_crit = Pt_in * (2 / (γ + 1))^(γ / (γ - 1))

    if p < p_crit
        # Choked flow
        mdot = Cd * A * Pt_in * sqrt(γ / (R * Tt_in)) *
               (2 / (γ + 1))^((γ + 1) / (2 * (γ - 1)))
    elseif p < Pt_in
        # Subsonic flow
        if !IsGas
            mdot = Cd * A * sqrt(2 * ρt * (Pt_in - p))
            println("p < Pt_in")
            println("p = $p, p_crit = $p_crit, Pt_in = $Pt_in")
        else
            ratio = p / Pt_in
            mdot = Cd * A * Pt_in * sqrt(γ / (R * Tt_in)) *
                   ratio^(1 / γ) *
                   sqrt((2 / (γ - 1)) * (1 - ratio^((γ - 1) / γ)))
        end
    else
        # Tank pressure has reached inlet pressure
        mdot *= 1 - Int(IsCoff)
    end
    mdot_vec[i] = mdot

    # Energy update
    U = m * R * T / (γ - 1)
    U += mdot * dt * γ * R * Tt_in / (γ - 1)

    # Update state variables
    m += mdot * dt
    T = (γ - 1) * U / (m * R)
    p = m * R * T / V

    # Check for crossover
    if !reached && mdot_vec[i] <= 0
        reached = true
        p_reach = p
        t_reach = time[i]
        mfr_reach = mdot_vec[i]
        mfr_0 = mdot_vec[1]
        v_0 = mfr_0 / (ρs * A360)

        println("Initial Inlet mass flow rate = $(mfr_0) kg/s at t = 0 s")
        println("Initial Inlet velocity = $(v_0) m/s at t = 0 s")
        println("-------------------------")
        println("Tank pressure p = $(round(p_reach)) Pa reached Pt_in = $(round(Pt_in)) Pa at t = $(t_reach) s")
        println("Final Inlet mass flow rate = $(mdot_vec[i]) kg/s at t = $(t_reach) s")
        # Optionally: break
    end
end

# Plotting
plot_layout = @layout [a; b; c; d]
p1 = plot(time, P ./ 1e5, lw=2, color=:blue, ylabel="Pressure [bar]", title="Tank Pressure Over Time", grid=true)
p2 = plot(time, Tvec, lw=2, color=:red, ylabel="Temperature [K]", title="Tank Temperature Over Time", grid=true)
p3 = plot(time, M, lw=2, color=:green, ylabel="Mass [kg]", title="Tank Mass Over Time", grid=true)
p4 = plot(time, mdot_vec, lw=2, color=:magenta, ylabel="Mass Flow Rate [kg/s]", xlabel="Time [s]", title="Mass Flow Rate Over Time", grid=true)
plot(p1, p2, p3, p4, layout=plot_layout, size=(800, 900))