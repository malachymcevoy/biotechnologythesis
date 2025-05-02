# Biotechnology Thesis

This GitHub contains code relevant to B223626's Biotechnology Honour's Project, entitled "Cell-ebrating Efficiency: Optimisation of Biological Systems".

Some files use Pluto, an interactive notebook environment for Julia (similar to Jupyter). To run these files, type the following into a Julia terminal:
using Pkg; Pkg.add(Pluto)
import Pluto; Pluto.run()

NOTE: Julia's ODE solvers are updated regularly, and the tolerances of solvers like RadauIIA5 are altered (often made stricter). This can lead to errors in running some of the code in this repository. Decreasing the abstol and reltol of the solvers should help with this issue. If this doesn't work, alternate solvers like Rodas5(), Rosenbrock23() and RadauIAA9() are all valid alternatives. For exact package specifics for when all this code was run in working order, see Project.toml and Manifest.toml.