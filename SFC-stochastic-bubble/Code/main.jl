using DifferentialEquations
using Plots
using Statistics
using CSV
using DataFrames
using Printf
using Random
using LaTeXStrings
using Printf
using Statistics, StatsBase
using FilePathsBase  # for path join if needed

#using CairoMakie

# Define the path where you want to save the file
#save_path_simu = raw"C:\Users\Adrien\Dropbox\Adrien & Matheus\asset price bubbles\experiments\CapeTown\simulation_results.png"
#save_path_par = raw"C:\Users\Adrien\Dropbox\Adrien & Matheus\asset price bubbles\experiments\CapeTown\parameters.csv"
#save_path_val = raw"C:\Users\Adrien\Dropbox\Adrien & Matheus\asset price bubbles\experiments\CapeTown\final_values.csv"
save_dir = raw"C:\Users\Adrien\Dropbox\Adrien & Matheus\asset price bubbles\experiments\CapeTown"

# Utility: Export parameters to CSV
function export_parameters_to_csv(par, filename)
    df = DataFrame(Parameter=String[], Value=Float64[])
    for (key, value) in pairs(par)
        push!(df, (string(key), value))
    end
    # Add long-term mean of μ
    mu_0 = par.r_L - 0.5 * par.σ^2
    push!(df, ("μ₀ (long-term mean)", mu_0))
    var_0 = par.σ^2/(2*par.η_μ)
    push!(df, ("Var₀ (long-term variance)", var_0))
    CSV.write(filename, df)
end



# Drift function
drift!(du, u, p, t) = begin
    ω, e, m, ℓ, S, μ= u # 7 dimensions
    
    κ_L      = p.r_L                                                # rate of reimbursement = interest rate
    r        = min(p.r_max, p.r_L + p.ρ_1*(exp(-p.ρ_2*(μ-p.r_L))))  # stochastic rate
    π_f      = 1 - ω  - p.δ*p.ν + p.r_M*m - r*ℓ                     # profit share
    Δ        = clamp(p.Δ_0 + p.Δ_1 * π_f, p.Δ_min, p.Δ_max)         # dividend function for distributed profits
    ι        = p.ηₚ * (p.ξ * ω - 1)                                  # inflation
    inv      = clamp(p.κ₀ + p.κ₁ * π_f,p.κ_min, p.κ_max)            # invested share of GDP
    Φ        = p.ϕ₀ + p.ϕ₁*e                                        # linear Phillips curve             
    g        = inv / p.ν - p.δ                                      # real growth rate
    nom_g    = g + ι                                                # nominal growth rate        
    Ψ        = clamp(p.ψ₁ * nom_g + p.ψ₀, p.ψ_min, p.ψ_max)         # Speculative function based on Ψ(g)
    Λ_p      = p.λ_p * max(0, Ψ)                                    # positive bubble trend from excess demand
    Λ_m      = p.λ_m * max(0, -Ψ)
    trend_p  = (p.J_p + log(1-p.J_p)) * Λ_p * p.η_μ - log(1-p.J_p)* Λ_p
    trend_m  = (log(1+p.J_m) - p.J_m) * Λ_m * p.η_μ - log(1+p.J_m)* Λ_m

    # Differential equations for (ω, e, m, l, S, μ)
    du[1] = ω * (Φ- p.α - (1 - p.γ) * ι)                               #ω
    du[2] = e * (g - p.α - p.β)                                        #e
    du[3] = π_f - (1-p.ζ)*(p.ν * g) + (r - κ_L)*ℓ - Δ + Ψ - nom_g * m  #m
    du[4] = p.ζ * (p.ν * g) + (r - κ_L)*ℓ + Ψ - nom_g*ℓ                #ℓ                                                      #f
    du[5] = S * (p.r_L + p.J_p * Λ_p - p.J_m * Λ_m)                    #S
    du[6] = p.η_μ * (p.r_L - 0.5*p.σ^2  - μ) + trend_p + trend_m       #μ
end

# Diffusion checked
diffusion!(du, u, p, t) = (du[5] = p.σ * u[5] ; du[6] = p.σ)

# rate of trend up for jumps down
function rate_down(u, p, t)
    ω, e, m, ℓ, S, μ = u

    r       = min(p.r_max, p.r_L + p.ρ_1 * exp(-p.ρ_2 * (μ - p.r_L)))
    π_f     = 1 - ω - p.δ * p.ν + p.r_M * m - r * ℓ
    ι       = p.ηₚ * (p.ξ * ω - 1)
    inv     = min(max(0, p.κ₀ + p.κ₁ * π_f), p.κ_max)
    g       = inv / p.ν - p.δ
    nom_g   = g + ι
    Ψ       = clamp(p.ψ₁ * nom_g + p.ψ₀, p.ψ_min, p.ψ_max)
    Λ_p     = p.λ_p * max(0, Ψ)
    return Λ_p
end

# rate for trend down to jump ups
function rate_up(u, p, t)
    ω, e, m, ℓ, S, μ = u

    r       = min(p.r_max, p.r_L + p.ρ_1 * exp(-p.ρ_2 * (μ - p.r_L)))
    π_f     = 1 - ω - p.δ * p.ν + p.r_M * m - r * ℓ
    ι       = p.ηₚ * (p.ξ * ω - 1)
    inv     = min(max(0, p.κ₀ + p.κ₁ * π_f), p.κ_max)
    g       = inv / p.ν - p.δ
    nom_g   = g + ι
    Ψ       = clamp(p.ψ₁ * nom_g + p.ψ₀, p.ψ_min, p.ψ_max)
    Λ_m     = p.λ_m * max(0, -Ψ)
    return Λ_m
end

function affect_down!(integrator)
    p = integrator.p
    integrator.u[5] *= (1.0 - p.J_p)
    integrator.u[6] += log(1 - p.J_p)
end

function affect_up!(integrator)
    p = integrator.p
    integrator.u[5] *= (1.0 + p.J_m)
    integrator.u[6] += log(1 + p.J_m)
end

# Problem builder with numerical choice of variable/const intensity of jump
function build_jump_problem(u0, tspan, par; jump_type=:variable)
    sde_prob = SDEProblem(drift!, diffusion!, u0, tspan, par, 
                          isoutofdomain=(u, p, t) -> any(x -> !isfinite(x), u))

    if jump_type == :variable
        jump_down = VariableRateJump(rate_down, affect_down!)
        jump_up = VariableRateJump(rate_up, affect_up!)
    else
        jump_down = ConstantRateJump((u, p, t) -> rate_down(u, p, t), affect_down!)
        jump_up = ConstantRateJump((u, p, t) -> rate_up(u, p, t), affect_up!)
    end

    return JumpProblem(sde_prob, Direct(), jump_down, jump_up)
end


# Run a single simulation
function run_simulation(par, tspan, u0, crisis_threshold, time_limit)
    prob = build_jump_problem(u0, tspan, par, jump_type=:constant)
    sol = nothing
    try
        sol = solve(prob, SRIW1(); isoutofdomain=(u, p, t) -> any(x -> !isfinite(x), u), maxiters=1e7)
    catch e1
        @warn "SRIW1 failed: $e1. Retrying with SOSRI."
        try
            sol = solve(prob, SOSRI(); isoutofdomain=(u, p, t) -> any(x -> !isfinite(x), u), maxiters=1e7)
        catch e2
            @warn "SOSRI also failed: $e2"
            return true, true  # count as both crisis and blowup
        end
    end
    is_blowup = sol.retcode != :Success || sol.t[end] < time_limit
    e_values = sol[2, :]
    m_values = sol[3, :]
    ℓ_values = sol[4, :]

    leverage_diffs = ℓ_values .- m_values  # = ℓ - m

    is_crisis = any(e -> e < crisis_threshold, e_values) ||
                any(df -> df > 10, leverage_diffs)
    #is_crisis = any(e -> e < crisis_threshold, e_values)
    return is_crisis, is_blowup
end


function monte_carlo(param_name::Symbol, param_values::Vector, par, num_sim, crisis_threshold, time_limit, tspan, u0)
    results = DataFrame(value=Float64[], prob_crisis=Float64[], prob_blowup=Float64[], prob_union=Float64[])
    for val in param_values
        par_updated = merge(par, NamedTuple{(param_name,)}((val,)))
        crisis_count = Threads.Atomic{Int}(0)
        blowup_count = Threads.Atomic{Int}(0)
        both_count = Threads.Atomic{Int}(0)

        @info "Running simulations for $param_name = $val"
        Threads.@threads for i in 1:num_sim
            crisis, blowup = run_simulation(par_updated, tspan, u0, crisis_threshold, time_limit)
            Threads.atomic_add!(crisis_count, crisis ? 1 : 0)
            Threads.atomic_add!(blowup_count, blowup ? 1 : 0)
            Threads.atomic_add!(both_count, (crisis && blowup) ? 1 : 0)
        end

        prob_crisis = crisis_count[] / num_sim
        prob_blowup = blowup_count[] / num_sim
        prob_both = both_count[] / num_sim
        prob_union = prob_crisis + prob_blowup - prob_both

        push!(results, (val, prob_crisis, prob_blowup, prob_union))
    end
    return results
end

#to simulate a single trajectory for plotting
function simulate_trajectory(par, tspan, u0; jump_type=:variable)
    prob = build_jump_problem(u0, tspan, par, jump_type=jump_type)
    sol = solve(prob, SRIW1())
    return sol
end


function plot_trajectory_and_save(par, tspan, u0, save_dir::String, basename::String)

    # Create save paths
    save_path_simu = joinpath(save_dir, basename * "_plot.png")
    save_path_par  = joinpath(save_dir, basename * "_params.csv")
    save_path_val  = joinpath(save_dir, basename * "_final_values.csv")

    r_from_mu(μ, p) = min(p.r_max, p.r_L + p.ρ_1 * exp(-p.ρ_2 * (μ - p.r_L)))
    sol = simulate_trajectory(par, tspan, u0)

    if sol === nothing || sol.retcode != :Success
        @warn "Simulation failed or returned invalid result."
        return
    end

    μ_vals = sol[6, :]
    r_values = [r_from_mu(μ, par) for μ in μ_vals]
    t_vals = sol.t

    # Compute discounted price: S̃ = S * exp(-r_L * t)
    S_vals = sol[5, :]
    S_discounted = [S * exp(-par.r_L * t) for (S, t) in zip(S_vals, t_vals)]

    f_values = let
        ω_vals, e_vals, m_vals, ℓ_vals = sol[1, :], sol[2, :], sol[3, :], sol[4, :]
        [begin
            r = r_from_mu(μ, par)
            π_f = 1 - ω - par.δ * par.ν + par.r_M * m - r * ℓ
            ι = par.ηₚ * (par.ξ * ω - 1)
            inv = min(max(0, par.κ₀ + par.κ₁ * π_f), par.κ_max)
            g = inv / par.ν - par.δ
            nom_g = g + ι
            clamp(par.ψ₁ * nom_g + par.ψ₀, par.ψ_min, par.ψ_max)
        end for (ω, e, m, ℓ, μ) in zip(ω_vals, e_vals, m_vals, ℓ_vals, μ_vals)]
    end

    net_charges = [r_values[i] * sol[4, i] - par.r_M * sol[3, i] for i in 1:length(r_values)]

    colors = ["#003f5c", "#7a5195", "#ef5675", "#ffa600"]

    p1 = Plots.plot(sol.t, sol[1, :], label = L"\omega", lw=2, xlabel = "Time", ylabel = "Level", title = L"\omega, e", color=colors[1])
    Plots.plot!(p1, sol.t, sol[2, :], label = L"e", lw=2, color=colors[3])

    p2 = Plots.plot(sol.t, net_charges, label = "net charges", lw=2, xlabel = "Time", ylabel = "Net Charges", title = "Net Charges and f", color=colors[1])
    p5 = twinx()
    Plots.plot!(p5, sol.t, f_values, label = L"f", lw=2, color=colors[4])

    # Modified plot with S̃ instead of S
    p3 = Plots.plot(sol.t, S_discounted, label = L"\tilde{S}", lw=2, xlabel = "Time", ylabel = L"\tilde{S}", title = L"\tilde{S}, \mu, r", color=colors[3])
    p4 = twinx()
    Plots.plot!(p4, sol.t, μ_vals, label = L"\mu", lw=2, ylabel = L"\mu,\, r", color=colors[4])
    Plots.plot!(p4, sol.t, r_values, label = L"r", lw=2, color=colors[1])

    plot_combined = plot(p1, p2, p3, layout = (3, 1), size=(900, 900))
    savefig(plot_combined, save_path_simu)
    export_parameters_to_csv(par, save_path_par)
    display(plot_combined)

    # Final values
    ω, e, m, ℓ, S, μ = sol.u[end]
    r = r_from_mu(μ, par)
    π_f = 1 - ω - par.δ * par.ν + par.r_M * m - r * ℓ
    ι = par.ηₚ * (par.ξ * ω - 1)
    inv = min(max(0, par.κ₀ + par.κ₁ * π_f), par.κ_max)
    g = inv / par.ν - par.δ
    nom_g = g + ι
    f = clamp(par.ψ₁ * nom_g + par.ψ₀, par.ψ_min, par.ψ_max)

    μ_mean = mean(μ_vals)
    μ_var = var(μ_vals)
    μ_skew = skewness(μ_vals)
    μ_kurt = kurtosis(μ_vals)

    μ_p = par.r_L - 0.5 * par.σ^2 + (log(1 - par.J_p) + par.J_p) * par.λ_p * f
    μ_m = par.r_L - 0.5 * par.σ^2 + (log(1 + par.J_m) - par.J_p) * par.λ_m * f

    variable_names = ["ω", "e", "m", "ℓ", "S", "μ", "f", "μ_mean", "μ_var", "μ_skew", "μ_kurt", "μ_p", "μ_m"]
    final_values = [ω, e, m, ℓ, S, μ, f, μ_mean, μ_var, μ_skew, μ_kurt, μ_p, μ_m]
    df_final = DataFrame(Variable = variable_names, Value = final_values)
    CSV.write(save_path_val, df_final)
end




function sweep_param_mc_plot(param_name::Symbol,
                             param_values::AbstractVector,
                             par,
                             num_simulations::Int,
                             crisis_threshold,
                             time_limit,
                             tspan,
                             u0;
                             xlabel = nothing,
                             ylabel::AbstractString = "Probability",
                             use_union::Bool = true,
                             xtick_mode::Symbol = :auto,   # :auto | :percent | :raw
                             xtick_every::Int = 3,
                             digits::Int = 4,
                             outdir::AbstractString = ".",
                             basename = nothing,
                             save_csv::Bool = true,
                             save_png::Bool = true,
                             show_plot::Bool = true)

    # 1) grid (arrondi optionnel)
    vals = round.(collect(param_values), digits = digits)

    # 2) Monte Carlo
    results = monte_carlo(param_name, vals, par, num_simulations,
                          crisis_threshold, time_limit, tspan, u0)  # :contentReference[oaicite:1]{index=1}

    # 3) Noms de fichiers
    base = basename === nothing ? "sweep_$(String(param_name))" : String(basename)
    csv_path = joinpath(outdir, "$(base).csv")
    png_path = joinpath(outdir, "$(base).png")

    # 4) Sauvegarde CSV
    if save_csv
        CSV.write(csv_path, results)
    end

    # 5) Choix de la série à tracer
    y = use_union ? results.prob_union : results.prob_crisis

    # 6) xticks
    xticks_tuple = nothing
    if xtick_mode != :auto
        idx = 1:xtick_every:length(vals)
        xticks_vals = vals[idx]
        if xtick_mode == :percent
            xticks_labels = [@sprintf("%.1f%%", 100 * x) for x in xticks_vals]
        elseif xtick_mode == :raw
            xticks_labels = string.(xticks_vals)
        else
            error("xtick_mode must be :auto, :percent, or :raw")
        end
        xticks_tuple = (xticks_vals, xticks_labels)
    end

    # 7) Labels
    xlab = xlabel === nothing ? String(param_name) : xlabel

    # 8) Plot
    p = plot(results.value, y;
             lw = 2,
             xlabel = xlab,
             ylabel = ylabel,
             title = "",
             label = false,
             grid = false,
             legend = false,
             framestyle = :box,
             xticks = xticks_tuple)

    if save_png
        savefig(p, png_path)
    end
    if show_plot
        display(p)
    end

    return results, p, (csv_path = csv_path, png_path = png_path)
end


# fonctionne avec des clés dynamiques (:r_L, :σ, etc.)
update_par(par, pairs::Pair{Symbol}...) = (; par..., pairs...)

function prob_union_at(par, num_simulations, crisis_threshold, time_limit, tspan, u0)
    crisis_count = Threads.Atomic{Int}(0)
    blowup_count = Threads.Atomic{Int}(0)
    both_count   = Threads.Atomic{Int}(0)

    Threads.@threads for i in 1:num_simulations
        crisis, blowup = run_simulation(par, tspan, u0, crisis_threshold, time_limit)
        Threads.atomic_add!(crisis_count, crisis ? 1 : 0)
        Threads.atomic_add!(blowup_count, blowup ? 1 : 0)
        Threads.atomic_add!(both_count, (crisis && blowup) ? 1 : 0)
    end

    pc = crisis_count[] / num_simulations
    pb = blowup_count[] / num_simulations
    p_union = pc + pb - both_count[] / num_simulations
    return p_union
end



function heatmap_prob_union(param1::Symbol, values1::AbstractVector,
                            param2::Symbol, values2::AbstractVector,
                            base_par, num_simulations, crisis_threshold,
                            time_limit, tspan, u0;
                            savepath::Union{String,Nothing}=nothing)

    n1, n2 = length(values1), length(values2)
    Z = Matrix{Float64}(undef, n1, n2)

    for (i, v1) in enumerate(values1)
        for (j, v2) in enumerate(values2)
            par_new = update_par(base_par, param1 => v1, param2 => v2)
            Z[i, j] = prob_union_at(par_new, num_simulations, crisis_threshold, time_limit, tspan, u0)
        end
    end

    p = heatmap(values2, values1, Z;
        xlabel = string(param2),
        ylabel = string(param1),
        title  = "P(crisis ∪ blowup)",
        colorbar_title = "Probability",
        c = cgrad([:lightyellow, :orange, :red, :black]),
        framestyle = :box
    )

    if savepath !== nothing
        savefig(p, savepath)
    end
    display(p)
    return Z, p
end






################################################################
################################################################

# Updated parameters for the system, including σ (sigma) for noise

par = (
    α = 0.02,          # alpha (1)
    β = 0.02,           # beta (2)
    γ = 0.9,            # money illusion parameter (3)
    δ = 0.04,           # delta (4)
    ν = 2.7,            # nu (5)
    r_L = 0.02,         # r_L (6) 
    r_M = 0.01,         # r_b for deposits (7)
    ηₚ = 0.192,         # eta_p (8)
    ξ = 1.875,          # xi (9)
    ϕ₀ = -0.292,        # phi_0 (10)
    ϕ₁ = 0.469,         # phi_1 (11)
    κ₀ = 0.0318,        # kappa_0 (12)
    κ₁ = 0.575,         # kappa_1 (13)
    κ_min = 0,          # kappa_min
    κ_max = 0.30,       # kappa_2 (14)
    ψ₀ = -0.075,        # psi_0 speculative level at nom_g=0
    ψ₁ = 3.75,          # psi_1 is slope of speculation wrt to nom_g
    ψ_min = -0.15,       # lower bound for speculation
    ψ_max = 0.3,        # upper bound for speculation UPDATED
    Δ_min = 0,          # min share of dividend
    Δ_max = 0.30,       # max share of dividend
    Δ_0 = -0.078,       # constant of affine share
    Δ_1 = 0.553,        # rate of affine share
    ζ = 0.8,            # proportion of investment paid with loans
    η_μ = 0.5,          # speed of adjustment of the market trend
    σ = 0.10,           # Volatility per year
    λ_p = 1.0,          # intensity rate for downward jumps
    λ_m = 1.0,         # intensity rate for upward jumps
    J_p = 0.1,         # downward jump size
    J_m = 0.1,            # upward jump size
    ρ_1 = 0.01,        # stochastic part of the interest rate
    ρ_2 = 5,           # stochasticpart of the interest rate exp
    κ_L = 0.02,         # rate of repayment on loans
    r_max = 0.2         # maximal interest rate
)
# Common settings
u0 = [0.7, 0.9, 0.3, 0.5, 1, 0.025]



# Helper: create a new NamedTuple like `par` but with some fields changed
update_par(par; kwargs...) = (; par..., kwargs...)

runs = [
    (tag = "411a", changes = (r_L = 0.020, λ_p = 0.0, λ_m = 0.0, ρ_1 = 0.0), tspan = (0.0, 300.0)),
    (tag = "411b", changes = (r_L = 0.15, λ_p = 0.0, λ_m = 0.0, ρ_1 = 0.0), tspan = (0.0, 100.0)),
    (tag = "411c", changes = (r_L = 0.020, λ_p = 1.0, λ_m = 1.0, ρ_1 = 0.0), tspan = (0.0, 300.0)),
    (tag = "411d", changes = (r_L = 0.15, λ_p = 1.0, λ_m = 1.0, ρ_1 = 0.0), tspan = (0.0, 100.0)),   
    (tag = "411e", changes = (r_L = 0.020, λ_p = 0.0, λ_m = 0.0, ρ_1 = 0.01), tspan = (0.0, 300.0)),
    (tag = "411f", changes = (r_L = 0.15, λ_p = 0.0, λ_m = 0.0, ρ_1 = 0.01), tspan = (0.0, 100.0)),
    (tag = "412a", changes = (), tspan = (0.0, 200.0)),
    (tag = "412b", changes = (σ = 0.25,), tspan = (0.0, 200.0)),
    (tag = "412c", changes = (η_μ = 5.0,), tspan = (0.0, 200.0)),
    (tag = "412d", changes = (η_μ = 0.2,), tspan = (0.0, 200.0)),   
    (tag = "412e", changes = (η_μ = 5.0, ρ_2 = 3), tspan = (0.0, 200.0)),
    (tag = "412f", changes = (η_μ = 0.2, ρ_2 = 3), tspan = (0.0, 200.0)),   
]

Random.seed!(666)
for r in runs
    par_run = update_par(par; r.changes...)
    plot_trajectory_and_save(par_run, r.tspan, u0, save_dir, r.tag)
end





# SETTING TO DEFINE CRISIS AND PROBA ESTIMATES
tspan = (0.0, 150.0)
num_simulations = 5000
crisis_threshold = 0.05
time_limit = 100.0

# Sweep r_L (ticks en % comme ton exemple)
r_L_values = range(0.01, stop = 0.05, length = 50)
results_rL, p_rL, paths_rL = sweep_param_mc_plot(:r_L, r_L_values, par,
    num_simulations, crisis_threshold, time_limit, tspan, u0;
    xlabel = L"r_L",
    xtick_mode = :percent,
    xtick_every = 3,
    digits = 4,
    outdir = save_dir,
    basename = "probability_vs_rL",
    use_union = true)

# Sweep η_μ (ticks “raw” comme ton exemple)
η_μ_values = range(0.03, stop = 0.7, length = 50)
results_eta, p_eta, paths_eta = sweep_param_mc_plot(:η_μ, η_μ_values, par,
    num_simulations, crisis_threshold, time_limit, tspan, u0;
    xlabel = L"\eta_{\mu}",
    xtick_mode = :raw,
    xtick_every = 5,
    digits = 3,
    outdir = save_dir,
    basename = "probability_vs_eta_mu",
    use_union = true)


# Sweep ρ_2 (ticks “raw” comme ton exemple)
ρ_2_values = range(3, stop = 15, length = 50)
results_rho2, p_rho2, paths_rho2 = sweep_param_mc_plot(:ρ_2, ρ_2_values, par,
    num_simulations, crisis_threshold, time_limit, tspan, u0;
    xlabel = L"\rho_{2}",
    xtick_mode = :raw,
    xtick_every = 5,
    digits = 3,
    outdir = save_dir,
    basename = "probability_vs_rho_2",
    use_union = true)



# Sweep σ (ticks “raw” comme ton exemple)
σ_values = range(0.10, stop = 0.30, length = 50)
results_sigma, p_sigma, paths_sigma = sweep_param_mc_plot(:σ, σ_values, par,
    num_simulations, crisis_threshold, time_limit, tspan, u0;
    xlabel = L"\sigma",
    xtick_mode = :raw,
    xtick_every = 5,
    digits = 3,
    outdir = save_dir,
    basename = "probability_vs_sigma",
    use_union = true)




tspan = (0.0, 120.0)
num_simulations = 600
crisis_threshold = 0.05
time_limit = 100.0

η_μ_vals= range(0.05, stop=10, length=70)
ρ_2_vals  = range(1, stop=25, length=70)

Z, p = heatmap_prob_union(:η_μ, η_μ_vals, :ρ_2, ρ_2_vals,
                          par, num_simulations, crisis_threshold,
                          time_limit, tspan, u0;
                          savepath = joinpath(save_dir, "heatmap_η_μ_rho_2.png"))


r_L_vals= range(0.01, stop=0.05, length=70)
σ_vals  = range(0.05, stop=0.25, length=70)

Z, p = heatmap_prob_union(:r_L, r_L_vals, :σ, σ_vals,
                          par, num_simulations, crisis_threshold,
                          time_limit, tspan, u0;
                          savepath = joinpath(save_dir, "heatmap_r_L_σ.png"))
