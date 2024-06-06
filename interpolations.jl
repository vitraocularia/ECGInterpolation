module ECGInterpolations

function estimate_derivatives(x, y)
    n = length(x)
    dydx = zeros(n)

    dydx[2:n-1] = (y[3:n] .- y[1:n-2]) ./ (x[3:n] .- x[1:n-2])

    dydx[1] = (y[2] - y[1]) / (x[2] - x[1])
    dydx[n] = (y[n] - y[n-1]) / (x[n] - x[n-1])

    return dydx
end

function select_nodes(x, y, num_nodes)
	n = length(x)
    step = n ÷ (num_nodes - 1)
    indices = [1:step:n; n]
    return x[indices], y[indices]
end

# HERMITE SPLINES INTERPOLATION

function hermite_basis(t)
	h00 = 2t^3 - 3t^2 + 1
    h10 = t^3 - 2t^2 + t
    h01 = -2t^3 + 3t^2
    h11 = t^3 - t^2
    return h00, h10, h01, h11
end

function hermite_spline_interpolation(x, y, n)

	x_nodes, y_nodes = select_nodes(x, y, n)
	dydx = estimate_derivatives(x_nodes, y_nodes)
	x_new = collect(range(minimum(x), stop=maximum(x), length=length(x)))
	y_new = zeros(length(x_new))

	for i in 1:(n-1)
		x0, x1 = x_nodes[i], x_nodes[i+1]
		y0, y1 = y_nodes[i], y_nodes[i+1]
		dy0, dy1 = dydx[i], dydx[i + 1]

		 for j in eachindex(x_new)
            xx = x_new[j]
            if x0 <= xx <= x1
                t = (xx - x0) / (x1 - x0)
                h00, h10, h01, h11 = hermite_basis(t)
                y_new[j] = h00 * y0 + h10 * (x1 - x0) * dy0 +
                           h01 * y1 + h11 * (x1 - x0) * dy1
            end
        end
	end
	return x_new, y_new
end

# B SPLINES INTERPOLATION

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

function b_spline_interpolation(x, y, num_nodes, k=3)
	
	x_nodes, y_nodes = select_nodes(x, y, num_nodes)

	n = length(x_nodes)
    knots = [x[1] for _ in 1:k]  # k repeated knots at the beginning
    append!(knots, x_nodes)
    append!(knots, [x_nodes[end] for _ in 1:k])  # k repeated knots at the end

	x_new = collect(range(minimum(x), stop=maximum(x), length=length(x)))
    y_new = zeros(length(x_new))

	for i in eachindex(x_new)
        for j in 1:n
            y_new[i] += y_nodes[j] * b_spline_basis(j, k, x_new[i], knots)
        end
    end
	
	return x_new, y_new
end

# CUBIC SPLINES INTERPOLATION

function cubic_spline_interpolation(x, y, num_nodes)
	t = collect(range(minimum(x), stop=maximum(x), length=num_nodes))
	n = length(x)
	h = diff(x)
	alpha = zeros(Float64, n-1)
	l = zeros(Float64, n)
	mu = zeros(Float64, n)
	z = zeros(Float64, n)
	c = zeros(Float64, n)
	b = zeros(Float64, n-1)
	d = zeros(Float64, n-1)
	
	for i in 2:n-1
		alpha[i] = (3/h[i])*(y[i+1]-y[i]) - (3/h[i-1])*(y[i]-y[i-1])
	end
	
	l[1] = 1.0
	mu[1] = 0.0
	z[1] = 0.0
	    
	for i in 2:n-1
		l[i] = 2*(x[i+1] - x[i-1]) - h[i-1]*mu[i-1]
		mu[i] = h[i]/l[i]
		z[i] = (alpha[i] - h[i-1]*z[i-1]) / l[i]
	end
	
	l[n] = 1.0
	z[n] = 0.0
	c[n] = 0.0
	
	for j in (n-1):-1:1
		c[j] = z[j] - mu[j]*c[j+1]
		b[j] = (y[j+1] - y[j])/h[j] - h[j]*(c[j+1] + 2*c[j])/3
		d[j] = (c[j+1] - c[j])/(3*h[j])
	end
	    
	function cubic_poly(j, t)
		return y[j] + b[j]*(t-x[j]) + c[j]*(t-x[j])^2 + d[j]*(t-x[j])^3
	end

	y_new = zeros(Float64, length(t))
	for i in eachindex(t)
		for j in 1:n-1
			if x[j] <= t[i] <= x[j+1]
				y_new[i] = cubic_poly(j, t[i])
			end
		end
	end
	
	return t, y_new
end

end