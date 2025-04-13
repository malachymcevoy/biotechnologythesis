### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 7e71dbe7-6493-42ea-9c27-ae3faf17e3b8
using PlutoUI, Plots, DifferentialEquations, Optim

# ╔═╡ bf9ef966-f4f1-11ef-0098-a1e60a9f658a
md"# Lotka-Volterra Equations"

# ╔═╡ 4caf3fef-0a72-4fed-936f-d86f925427c8
md"""
The Lotka–Volterra equations are a pair of first-order nonlinear differential equations used to describe the dynamics of biological systems when two species interact as predator and as prey. 

The populations change through time according to the pair of equations:

$\frac {dx}{dt} = \alpha x - \beta x y$
$\frac {dy}{dt} = -\gamma y + \delta x y$

Where the variables $x$ and $y$ are the population density of prey and predator respectively, $\alpha$ and $\beta$ describe the maximum prey per capita growth rate and the effect of the presence of predators on the prey death rate respectively, and $\gamma$ and $\delta$ describe respectively the predator's per capita death rate and the effect of the presence of prey on the predator's growth rate.
"""

# ╔═╡ 72c4fddb-fb8f-491f-9e93-91a07990a869
md"## Plotting the equations"

# ╔═╡ 393eb9a6-1f39-4dc8-bfbf-aa6a76c64568
md"Import Packages:"

# ╔═╡ 347d8e7d-72b3-485d-a1a8-c0f93ec16e83
md"Define the function:"

# ╔═╡ 17aa8d2f-85cb-4840-aa1d-78497179692c
function lotka_volterra(du, u, p, t)
    α, β, γ, δ = p
    x, y = u
    du[1] = α * x - β * x * y   # Prey population equation
    du[2] = -γ * y + δ * x * y  # Predator population equation
end

# ╔═╡ bc44fe5f-0873-49f6-96e1-f6a425c31fdd
md"Define the initial values and set the parameters"

# ╔═╡ cb96545b-cbb8-4e74-a73a-0bc826e9c5ae
begin
	#plt2 = plot(sol[1, :], sol[2, :], linewidth=2, xlabel="Prey Population", ylabel="Predator Population",
		            #title="Phase Plane", legend=false, size=(600,400))
end

# ╔═╡ 7b9e1a27-a0a8-4d5e-9a44-fb5c412da813
md"Adjust the sliders below to change the parameters and initial values:"

# ╔═╡ a186903a-c124-4c78-b094-831987b6c92b
@bind reset Button("Reset")

# ╔═╡ d2dfc322-dd34-46ab-9991-a76985f0acb1
begin
	reset
	α = @bind α Slider(0.0:0.01:2.0, default = 0.5, show_value=true)
	β = @bind β Slider(0.1:0.01:1.0, default = 1.0, show_value=true)
	γ = @bind γ Slider(0.1:0.01:3.0, default = 3.0, show_value=true)
	δ = @bind δ Slider(0.1:0.01:2.0, default = 1.0, show_value=true)
	prey_init = @bind prey_init Slider(1:1:50, default=5, show_value=true)
	predator_init = @bind predator_init Slider(1:1:20, default = 5, show_value=true)
	
	md"""
	- **α (Prey Birth Rate)**: $α  
	- **β (Predation Rate)**: $β  
	- **γ (Predator Death Rate)**: $γ  
	- **δ (Predator Birth Rate)**: $δ  
	- **Initial Prey Population**: $prey_init  
	- **Initial Predator Population**: $predator_init
	"""
end

# ╔═╡ bc78147b-50a9-465e-969a-4d0191d4f6eb
begin
	reset
	md"- **Years**: $(@bind t_end Slider(0:100, 100, true))"
end

# ╔═╡ 071fdad1-0c47-42e7-b97d-795f4c4aa641
begin
	init = [prey_init, predator_init]  
	params = [α, β, γ, 0.005093643076928691] 
	t_begin = 0
	tspan = (t_begin, t_end)
end

# ╔═╡ 448b6e32-e745-4f18-bba1-32fb9ece8b06
begin
	prob = ODEProblem(lotka_volterra, init, tspan, params, isoutofdomain=(u, p, t) -> any(x -> x < 0, u))
	sol = solve(prob, Tsit5())
	
	plt1 = plot(sol.t, sol[1, :], label="Prey", linewidth=2, xlabel="Time", ylabel="Population", title="Population vs. Time", size=(600,400), xlims=(0,t_end))
	plot!(sol.t, sol[2, :], label="Predator", linewidth=2)
	
	plot(plt1, size=(700,350))  
end

# ╔═╡ 44971249-97c7-4b66-90f8-7ee2d5322709
md"##### Return maximum prey and maximum predator values (relevant for later optimization)"

# ╔═╡ fd4aa34a-d5fa-4006-baee-4342f9922478
begin
	peak_prey = maximum(sol[1, :])
	println("Peak prey number = ", peak_prey)
end

# ╔═╡ 3bff9c0c-9cb1-4b96-8af7-bc4cecb0a9bb
begin
	peak_pred = maximum(sol[2, :])
	println("Peak predator number = ", peak_pred)
end

# ╔═╡ 100517ad-1b25-4560-b25b-5ead63198744
md"# Phase plot"

# ╔═╡ 2d1ca916-4712-47e2-b8a1-c3703248d7bc
phase = plot(sol, idxs = (1,2),
	legend = false,
	linewidth = 2,
	title = "Lotka_volterra Equations (Phase Space)",
	xaxis = "Prey Population",
	yaxis = "Predator Population",
	formatter = :plain,
	widen = true,
	xlims = (0.0, 30.0),
	ylims = (0.0, 15.0),
	aspect_ratio = 1.1
);

# ╔═╡ 84465733-c14a-442c-b98c-5e0f64df45e4
scatter!(phase, (sol(t_end)[1], sol(t_end)[2]), color = :red, markersize = 5, xlims=(0, 100), ylims=(0,50))

# ╔═╡ Cell order:
# ╟─bf9ef966-f4f1-11ef-0098-a1e60a9f658a
# ╟─4caf3fef-0a72-4fed-936f-d86f925427c8
# ╟─72c4fddb-fb8f-491f-9e93-91a07990a869
# ╟─393eb9a6-1f39-4dc8-bfbf-aa6a76c64568
# ╠═7e71dbe7-6493-42ea-9c27-ae3faf17e3b8
# ╟─347d8e7d-72b3-485d-a1a8-c0f93ec16e83
# ╠═17aa8d2f-85cb-4840-aa1d-78497179692c
# ╟─bc44fe5f-0873-49f6-96e1-f6a425c31fdd
# ╠═071fdad1-0c47-42e7-b97d-795f4c4aa641
# ╠═448b6e32-e745-4f18-bba1-32fb9ece8b06
# ╟─cb96545b-cbb8-4e74-a73a-0bc826e9c5ae
# ╟─7b9e1a27-a0a8-4d5e-9a44-fb5c412da813
# ╟─a186903a-c124-4c78-b094-831987b6c92b
# ╠═d2dfc322-dd34-46ab-9991-a76985f0acb1
# ╟─bc78147b-50a9-465e-969a-4d0191d4f6eb
# ╟─44971249-97c7-4b66-90f8-7ee2d5322709
# ╟─fd4aa34a-d5fa-4006-baee-4342f9922478
# ╟─3bff9c0c-9cb1-4b96-8af7-bc4cecb0a9bb
# ╟─100517ad-1b25-4560-b25b-5ead63198744
# ╟─2d1ca916-4712-47e2-b8a1-c3703248d7bc
# ╠═84465733-c14a-442c-b98c-5e0f64df45e4
