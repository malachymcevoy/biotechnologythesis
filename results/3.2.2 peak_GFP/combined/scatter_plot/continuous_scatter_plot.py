import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm

def model(t, u, p):
    qx, ng, ql, nG, ein, mu, Din, kb, ku, kd, km, Rin, ktx, Ktx, ktl, Ktl = p
    e, D, m, cl, G, Gm, R = u
    vtx = D * (ktx / ng) * (e / (Ktx + e))
    vtl = cl * (ktl / nG) * (e / (Ktl + e))
    du = np.zeros(7)
    du[0] = - (vtx * qx * ng + vtl * ql * nG) + ein - mu * e
    du[1] = Din - mu * D
    du[2] = vtx - kb * R * m + ku * cl + vtl - kd * m - mu * m
    du[3] = kb * R * m - ku * cl - vtl - kd * cl - mu * cl
    du[4] = vtl - km * G - mu * G
    du[5] = km * G - mu * Gm
    du[6] = -kb * R * m + ku * cl + vtl + kd * cl + Rin - mu * R
    return du

def run_simulation(dna_conc, energy_conc, tau):
    initial_conditions = [energy_conc, dna_conc, 0.0, 0.0, 0.0, 0.0, 1.51]
    p = [2.0, 833.0, 4.0, 236.0, 33600.0 * (0.2 / tau), (0.2 / tau), 0.005 * (0.2 / tau),
         1000.0, 1.0, 0.038153465, 0.084, 1.51 * (0.2/tau), 3750.0, 54.75, 105.0, 100.0]
    
    sol = solve_ivp(
        fun=lambda t, u: model(t, u, p),
        t_span=(0, 1000), 
        y0=initial_conditions,
        method='Radau',
        t_eval=np.linspace(0, 1000, 100)
    )
    
    return np.max(sol.y[5]) if sol.success else np.nan

def create_3d_scatter_plot():
    dna_vals = np.linspace(0.001, 0.05, 20)
    energy_vals = np.linspace(30000, 60000, 20)  
    tau_vals = np.linspace(10, 180, 20)  
    
    raw_points = []
    gfp_vals = []
    
    total_simulations = len(dna_vals) * len(energy_vals) * len(tau_vals)
    print(f"Running {total_simulations} simulations...")
    
    with tqdm(total=total_simulations) as pbar:
        for dna in dna_vals:
            for energy in energy_vals:
                for tau in tau_vals:
                    gfp = run_simulation(dna, energy, tau)
                    if not np.isnan(gfp):
                        raw_points.append((dna, energy, tau))
                        gfp_vals.append(gfp)
                    pbar.update(1)
    
    print(f"Generated {len(gfp_vals)} valid data points")
    
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')
    
    dna = [point[0] for point in raw_points]
    energy = [point[1] for point in raw_points]
    tau = [point[2] for point in raw_points]
    
    norm = plt.Normalize(vmin=0, vmax=np.max(gfp_vals))
    cmap = plt.cm.viridis
    
    sc = ax.scatter(dna, energy, tau, c=gfp_vals, cmap=cmap, norm=norm)
    
    ax.set_xlabel('DNA Concentration (µM)')
    ax.set_ylabel('Energy Concentration (µM)')
    ax.set_zlabel('Dilution Interval τ (min)')
    ax.set_title('3D Scatter Plot of GFP Production')
    
    cbar = plt.colorbar(sc, ax=ax, pad=0.1, shrink=0.7)
    cbar.set_label('Peak GFP Level')
    
    max_idx = np.argmax(gfp_vals) if gfp_vals else None
    if max_idx is not None:
        ax.text(dna[max_idx], energy[max_idx], tau[max_idx],
                f'Max GFP: {gfp_vals[max_idx]:.2f}',
                fontsize=10, color='red')
    
    plt.savefig('gfp_3d_scatter_plot.png', dpi=400)
    plt.show()

create_3d_scatter_plot()