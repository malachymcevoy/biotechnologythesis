# Biotechnology Thesis

This GitHub contains code relevant to my Biotechnology Honour's Project, entitled "Applying Optimization Algorithms to Biological Systems and Cell-Free Models".

Some files use Pluto, an interactive notebook environment for Julia (similar to Jupyter). To run these files, type the following into a Julia terminal:
using Pkg; Pkg.add(Pluto)
import Pluto; Pluto.run()

NOTE: Julia's ODE solvers are update regularly, and the tolerances of solvers like RadauIIA5 are altered (often made stricter). This can lead to errors in running some of the code in this repository. Decreasing the abstol and reltol of the solvers should help with this issue