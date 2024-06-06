### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ e06e122c-1f53-11ef-19ee-8b0eab47a233
begin
	using Pkg
	Pkg.add(["CSV", "DataFrames", "Plots", "PlutoUI"])
	using CSV, DataFrames, Plots, PlutoUI
end

# ╔═╡ 9b8f1876-81c7-43a4-b774-705a71d3a2aa
include("interpolations.jl")

# ╔═╡ 293dd522-bae7-463f-b7e1-7a43b69d49c6
function estimate_derivatives(x, y)
    n = length(x)
    dydx = zeros(n)

    dydx[2:n-1] = (y[3:n] .- y[1:n-2]) ./ (x[3:n] .- x[1:n-2])

    dydx[1] = (y[2] - y[1]) / (x[2] - x[1])
    dydx[n] = (y[n] - y[n-1]) / (x[n] - x[n-1])

    return dydx
end

# ╔═╡ 5f228002-b1e6-4afc-9a1a-0bff64d15622
function hermite_basis(t)
	h00 = 2t^3 - 3t^2 + 1
    h10 = t^3 - 2t^2 + t
    h01 = -2t^3 + 3t^2
    h11 = t^3 - t^2
    return h00, h10, h01, h11
end

# ╔═╡ b87408b1-b51c-46ac-8ed9-9f16eea84049


# ╔═╡ 046dea7c-3b1b-4b92-a7e6-f1bb5364b41b
begin
	# Preparing some data
	mitbih_test = CSV.File("./mitbih_test.csv", header=0) |> DataFrame
	
	# Usuwam ostatnią kolumnę, bo zawiera wynik klasyfikacji
	select!(mitbih_test, Not(size(mitbih_test)[2])) 
	data = Matrix(mitbih_test)
end

# ╔═╡ 3f8bd629-72e7-4583-a1af-0bc57db20bc2
md"n: $(@bind n PlutoUI.Slider(2:188; default=188))"

# ╔═╡ c6d30402-92e7-4148-8462-849c9774a14f
length(data[5,:])

# ╔═╡ c03f3095-7015-4a12-ad27-c530ceeced0c
# trying another type of splines - B-splines

# ╔═╡ 9759c536-9f4e-43a4-a27d-2a796adc4654
function b_spline_basis(j, k, t, nodes)
	if k == 1
		return nodes[j] <= t < nodes[j+1] ? 1.0 : 0.0
	else
		l = 0.0
		r = 0.0

		if (nodes[j+k-1] - nodes[j]) != 0
			l = (t - nodes[j]) / (nodes[j+k-1] - nodes[j]) * b_spline_basis(j, k-1, t, nodes)
		end

		if (nodes[j+k] - nodes[j+1]) != 0
			r = (nodes[j+k] - t) / (nodes[j+k] - nodes[j+1]) * b_spline_basis(j+1, k-1, t, nodes)
		end
		return l + r
	end
end

# ╔═╡ 90e31d88-6eb6-4007-9a35-0c093124e9e7
md"m: $(@bind m PlutoUI.Slider(2:188; default=188))"

# ╔═╡ ca3b7720-a31a-4fe3-9b2b-44bfb56b4c19


# ╔═╡ 068286d8-0548-42ad-90de-0afb5d8415ff
function select_nodes(x, y, num_nodes)
	n = length(x)
    step = n ÷ (num_nodes - 1)
    indices = [1:step:n; n]  # Ensure the last node is included
    return x[indices], y[indices]
end

# ╔═╡ 75bf1be3-711f-4f41-ab4c-8bc75e9fcb90
function hermite_spline_interpolation(x, y, n)

	x_nodes, y_nodes = select_nodes(x, y, n)
	dydx = estimate_derivatives(x_nodes, y_nodes)
	x_new = collect(range(minimum(x), stop=maximum(x), length=length(x)))
	y_new = zeros(length(x_new))

	for i in 1:(n-1)
		x0, x1 = x_nodes[i], x_nodes[i+1]
		y0, y1 = y_nodes[i], y_nodes[i+1]
		dy0, dy1 = dydx[i], dydx[i + 1]

		 for j in 1:length(x_new)
            xx = x_new[j]
            if x0 <= xx <= x1
                t = (xx - x0) / (x1 - x0)
                h00, h10, h01, h11 = hermite_basis(t)
                y_new[j] = h00 * y0 + h10 * (x1 - x0) * dy0 +
                           h01 * y1 + h11 * (x1 - x0) * dy1
            end
        end
	end
	return x_new, y_new, x_nodes, y_nodes
end

# ╔═╡ b5d3209f-1cfb-4fa0-b2c4-c04db5349b06
function plot_signal(amplitude, num_nodes)
	time = collect(1:length(amplitude))
	x_new, y_new, x_nodes, y_nodes = hermite_spline_interpolation(time, amplitude, num_nodes)
	plot(time, amplitude, label="Original Signal", lw=2, color=:black, alpha=0.5)
	plot!(x_new, y_new, label="Hermite Spline", lw=2, color=:red)
	plot!(title="Hermite Interpolation (n = $num_nodes)")
	scatter!(x_nodes, y_nodes, label="nodes", ms=2, color=:black)
end

# ╔═╡ b169183d-eea4-4fb3-a556-171dda0fb1f0
plot_signal(data[5,:], n)

# ╔═╡ 90a7dc70-7c30-42b7-b8bc-3e6970ce62ce
function b_spline_interpolation(x, y, num_nodes, k=3)
	
	x_nodes, y_nodes = select_nodes(x, y, num_nodes)

	n = length(x_nodes)
    knots = [x[1] for _ in 1:k]  # k repeated knots at the beginning
    append!(knots, x_nodes)
    append!(knots, [x_nodes[end] for _ in 1:k])  # k repeated knots at the end

	x_new = collect(range(minimum(x), stop=maximum(x), length=length(x)))
    y_new = zeros(length(x_new))

	for i in 1:length(x_new)
        for j in 1:n
            y_new[i] += y_nodes[j] * b_spline_basis(j, k, x_new[i], knots)
        end
    end
	
	return x_new, y_new, x_nodes, y_nodes
end

# ╔═╡ 0b3ad1ca-8e3f-456d-bc7f-0b05b52cc00d
function plot_signal2(amplitude, num_nodes)
	time = collect(1:length(amplitude))
	x_new, y_new, x_nodes, y_nodes = b_spline_interpolation(time, amplitude, num_nodes)
	plot(time, amplitude, label="Original Signal", lw=2, color=:black, alpha=0.5)
	plot!(x_new, y_new, label="B-Spline", lw=2, color=:red)
	plot!(title="B-Spline Interpolation (n = $num_nodes)")
	scatter!(x_nodes, y_nodes, label="nodes", ms=2, color=:black)
end

# ╔═╡ 90aaf451-4d44-4fa1-83c4-1bc116611dce
plot_signal2(data[5,:], m)

# ╔═╡ 49b3124b-7fae-4741-a694-2015435fcd1f
function cubic_spline_interpolation(x, y, n)
	t = collect(range(minimum(x), stop=maximum(x), length=n))
	n = length(x)
	h = diff(x)
	α = zeros(Float64, n-1)
	l = zeros(Float64, n)
	μ = zeros(Float64, n)
	z = zeros(Float64, n)
	c = zeros(Float64, n)
	b = zeros(Float64, n-1)
	d = zeros(Float64, n-1)
	
	for i in 2:n-1
		α[i] = (3/h[i])*(y[i+1]-y[i]) - (3/h[i-1])*(y[i]-y[i-1])
	end
	
	l[1] = 1.0
	μ[1] = 0.0
	z[1] = 0.0
	    
	for i in 2:n-1
		l[i] = 2*(x[i+1] - x[i-1]) - h[i-1]*μ[i-1]
		μ[i] = h[i]/l[i]
		z[i] = (α[i] - h[i-1]*z[i-1]) / l[i]
	end
	
	l[n] = 1.0
	z[n] = 0.0
	c[n] = 0.0
	
	for j in (n-1):-1:1
		c[j] = z[j] - μ[j]*c[j+1]
		b[j] = (y[j+1] - y[j])/h[j] - h[j]*(c[j+1] + 2*c[j])/3
		d[j] = (c[j+1] - c[j])/(3*h[j])
	end
	    
	function cubic_poly(j, t)
		return y[j] + b[j]*(t-x[j]) + c[j]*(t-x[j])^2 + d[j]*(t-x[j])^3
	end

	y_new = zeros(Float64, length(t))
	for i in 1:length(t)
		for j in 1:n-1
			if x[j] <= t[i] <= x[j+1]
				y_new[i] = cubic_poly(j, t[i])
			end
		end
	end
	
	return t, y_new
end

# ╔═╡ 5e047e5d-e1a0-4621-ba43-11e5384dcb0a
# ╠═╡ disabled = true
#=╠═╡
function evaluate_cubic_spline(x, x_nodes, b, c, d)
    n = length(x_nodes)
    y = similar(x)
    
    for i in 1:length(x)
        idx = searchsortedlast(x_nodes, x[i])
        if idx == 0
            idx = 1
        elseif idx == n
            idx -= 1
        end
        h = x[i] - x_nodes[idx]
        y[i] = b[idx] + (c[idx] * h) + (d[idx] * h^2)
    end
    return y
end

  ╠═╡ =#

# ╔═╡ 03bb687e-ac64-414a-84db-3772d866a4fe
function plot_cubic(amplitude, num_nodes)
	time = collect(1:length(amplitude))
	x_new, y_new = cubic_spline_interpolation(time, amplitude, num_nodes)
	plot(time, amplitude, label="Original Signal", lw=2, color=:black, alpha=0.4)
	plot!(x_new, y_new, label="Cubic Spline", lw=2, color=:red)
	plot!(title="Cubic Spline Interpolation")
end

# ╔═╡ 5e4c8b67-1fad-49e4-9307-870b82c0c745
md"o: $(@bind o PlutoUI.Slider(2:188; default=188))"

# ╔═╡ a1f205ea-8f5c-49c8-94d7-9d73d2a79d64
plot_cubic(data[5,:], o)

# ╔═╡ 2db57740-4b27-49b9-86ae-3141ac4d41f8


# ╔═╡ Cell order:
# ╠═e06e122c-1f53-11ef-19ee-8b0eab47a233
# ╠═293dd522-bae7-463f-b7e1-7a43b69d49c6
# ╠═5f228002-b1e6-4afc-9a1a-0bff64d15622
# ╠═b87408b1-b51c-46ac-8ed9-9f16eea84049
# ╠═75bf1be3-711f-4f41-ab4c-8bc75e9fcb90
# ╠═046dea7c-3b1b-4b92-a7e6-f1bb5364b41b
# ╠═b5d3209f-1cfb-4fa0-b2c4-c04db5349b06
# ╠═3f8bd629-72e7-4583-a1af-0bc57db20bc2
# ╠═b169183d-eea4-4fb3-a556-171dda0fb1f0
# ╠═c6d30402-92e7-4148-8462-849c9774a14f
# ╠═c03f3095-7015-4a12-ad27-c530ceeced0c
# ╠═9759c536-9f4e-43a4-a27d-2a796adc4654
# ╠═90a7dc70-7c30-42b7-b8bc-3e6970ce62ce
# ╠═90e31d88-6eb6-4007-9a35-0c093124e9e7
# ╠═0b3ad1ca-8e3f-456d-bc7f-0b05b52cc00d
# ╠═90aaf451-4d44-4fa1-83c4-1bc116611dce
# ╠═ca3b7720-a31a-4fe3-9b2b-44bfb56b4c19
# ╠═068286d8-0548-42ad-90de-0afb5d8415ff
# ╠═49b3124b-7fae-4741-a694-2015435fcd1f
# ╠═5e047e5d-e1a0-4621-ba43-11e5384dcb0a
# ╠═03bb687e-ac64-414a-84db-3772d866a4fe
# ╠═a1f205ea-8f5c-49c8-94d7-9d73d2a79d64
# ╠═5e4c8b67-1fad-49e4-9307-870b82c0c745
# ╠═9b8f1876-81c7-43a4-b774-705a71d3a2aa
# ╠═2db57740-4b27-49b9-86ae-3141ac4d41f8
