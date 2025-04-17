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
    params = [2.0, 833.0, 4.0, 236.0, energy_conc * 0.0, 0.0, dna_conc * 0.0,
              1000.0, 1.0, 0.038153465, 0.084, 1.51 * 0.0, 3750.0, 54.75, 105.0, 100.0]

    idx_e, idx_D_G, idx_m_G, idx_cl_G, idx_G, idx_G_m, idx_R = 1, 2, 3, 4, 5, 6, 7
    
    INDEX_DISC_A = [idx_e, idx_D_G, idx_R]
    CONC_DISC_A = [energy_conc, dna_conc, 1.51]
    INDEX_DISC_B = [idx_e, idx_R]
    CONC_DISC_B = [energy_conc, 1.51]
    INDEX_DISC_C = [idx_e, idx_R]
    CONC_DISC_C = [energy_conc, 1.51]

    dilSpecs = (INDEX_DISC_A, CONC_DISC_A, INDEX_DISC_B, CONC_DISC_B, INDEX_DISC_C, CONC_DISC_C)
    f_V = 0.2
    s1, s2, tmax, tsave, nspecies = 9.0, 0.0, 20.0, 100.0, 7.0

    function discreteDilute!(integrator)
        f_V = integrator.p[end - 6]
        
        if integrator.t < s1 * 60.0
            apply_dilution!(integrator, INDEX_DISC_A, CONC_DISC_A, f_V)
        elseif s1 * 60.0 <= integrator.t < s2 * 60.0
            apply_dilution!(integrator, INDEX_DISC_B, CONC_DISC_B, f_V)
        else
            apply_dilution!(integrator, INDEX_DISC_C, CONC_DISC_C, f_V)
        end
    end
    
    function apply_dilution!(integrator, indices, concentrations, f_V)
        for j in 1:7
            integrator.u[j] *= (1.0 - f_V)
        end
        for j in 1:length(indices)
            integrator.u[indices[j]] += f_V * concentrations[j]
        end
    end

    function get_callback(tau)
        return PeriodicCallback(discreteDilute!, tau; save_positions=(false, false))
    end

    params_flat = vcat(params, [f_V, tau, s1, s2, tmax, tsave, nspecies])
    periodcb = get_callback(tau)

    prob = ODEProblem(model!, initial_conditions, (0.0, 1000.0), params_flat)
    sol = solve(prob, RadauIIA5(), callback=periodcb, saveat=5.0)

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

    if !isempty(gfp_vals)
        max_idx = argmax(gfp_vals)
        max_point = raw_points[max_idx]
        println("Max GFP: $(round(gfp_vals[max_idx], digits=4))")
        println("DNA Concentration: $(round(max_point[1], digits=4)) µM")
        println("Energy Concentration: $(round(max_point[2], digits=4)) µM")
        println("Dilution Interval τ: $(round(max_point[3], digits=4)) min")
    end

    fig = scatter(
        [p[1] for p in raw_points], 
        [p[2] for p in raw_points], 
        [p[3] for p in raw_points],
        zcolor=gfp_vals,
        xlabel="DNA Concentration (µM)", 
        ylabel="                                    Energy Concentration (µM)", 
        # Add spaces to move label right
        zlabel="Dilution Interval τ (min)",
        title="3D Scatter Plot of GFP Production",
        markersize=5,
        markerstrokecolor=:auto,
        markerstrokewidth=0.2,
        c=:viridis,
        size=(1200, 900),
        label=false
    )

    savefig(fig, "gfp_3d_scatter_plot.png")
    display(fig)
end

create_3d_scatter_plot()
