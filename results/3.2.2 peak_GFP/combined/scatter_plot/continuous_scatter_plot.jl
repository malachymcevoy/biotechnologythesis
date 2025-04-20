using DifferentialEquations
using LinearAlgebra
using Plots
using ProgressMeter

function model!(du, u, p, t)
    qx, ng, ql, nG, ein, mu, Din, kb, ku, kd, km, Rin, ktx, Ktx, ktl, Ktl = p
    e, D, m, cl, G, Gm, R = u

    vtx = D * (ktx / ng) * (e / (Ktx + e))
    vtl = cl * (ktl / nG) * (e / (Ktl + e))

    du[1] = - (vtx * qx * ng + vtl * ql * nG) + ein - mu * e
    du[2] = Din - mu * D
    du[3] = vtx - kb * R * m + ku * cl + vtl - kd * m - mu * m
    du[4] = kb * R * m - ku * cl - vtl - kd * cl - mu * cl
    du[5] = vtl - km * G - mu * G
    du[6] = km * G - mu * Gm
    du[7] = -kb * R * m + ku * cl + vtl + kd * cl + Rin - mu * R
end

function run_simulation(dna_conc, energy_conc, tau)
    initial_conditions = [energy_conc, dna_conc, 0.0, 0.0, 0.0, 0.0, 1.51]
    p = [2.0, 833.0, 4.0, 236.0, 33600.0 * (0.2 / tau), (0.2 / tau), 0.005 * (0.2 / tau),
         1000.0, 1.0, 0.038153465, 0.084, 1.51 * (0.2/tau), 3750.0, 54.75, 105.0, 100.0]
    
    prob = ODEProblem(model!, initial_conditions, (0.0, 1000.0), p)
    sol = solve(prob, RadauIIA5(), saveat=100)

    return maximum(sol[6, :])
end

function create_3d_scatter_plot()
    dna_vals = range(0.001, 0.05, length=20)
    energy_vals = range(30000, 60000, length=20)
    tau_vals = range(10, 180, length=20)

    raw_points = []
    gfp_vals = Float64[]

    total_simulations = length(dna_vals) * length(energy_vals) * length(tau_vals)
    println("Running $total_simulations simulations...")

    @showprogress for dna in dna_vals, energy in energy_vals, tau in tau_vals
        gfp = run_simulation(dna, energy, tau)
        if !isnan(gfp)
            push!(raw_points, (dna, energy, tau))
            push!(gfp_vals, gfp)
        end
    end

    println("Generated $(length(gfp_vals)) valid data points")

    fig = scatter(
        [p[1] for p in raw_points], 
        [p[2] for p in raw_points], 
        [p[3] for p in raw_points],
        zcolor=gfp_vals,
        xlabel="DNA Concentration (µM)", ylabel="Energy Concentration (µM)",
        zlabel="Dilution Interval τ (min)",
        title="3D Scatter Plot of GFP Production",
        markersize=5,
        markerstrokecolor=:auto,
        markerstrokewidth=0.2,
        c=:viridis,
        size=(1200, 900)
    )

    if !isempty(gfp_vals)
        max_idx = argmax(gfp_vals)
        annotate!(fig, raw_points[max_idx][1], raw_points[max_idx][2], raw_points[max_idx][3],
                  text("Max GFP: $(round(gfp_vals[max_idx], digits=2))", 10, :red))
    end

    savefig(fig, "gfp_3d_scatter_plot.png")
    display(fig)
end

create_3d_scatter_plot()