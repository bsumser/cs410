using Plots
include("euler_methods.jl")

# u_t = κ*u_xx + F(t, x) on the domain 0 ≤ x ≤ 1, 0 ≤ t ≤ Tf
# u(0, x) = f(x)
# u(t, 0) = u(t, 1) = 0
# using second-order finite difference approx in space, and forward Euler in time.
# This becomes y' = A*y + b(t), where A is the discrete Laplacian matrix, and
# y0 = f(interior_spatial_nodes)

# λ = Δt/Δx^2, Stability → λ ≤ 1/(2κ)

κ = 2 # constant for thermal diffusivity
Δx = 0.1
λ = 0.5*(1/(2κ))
Δt = round(λ*Δx^2, digits = 10)

Tf = 2
M = Integer(Tf/Δt) # how many time steps to take

N = Integer(1/Δx) # N+1 total spatial nodes

#x = [0:Δx:1]
x = collect(range(0, 1, step = Δx))
x_int = x[2:end-1] #interior spatial nodes

t = collect(range(0, Tf, step = Δt))

A = zeros(N-1, N-1)
for i = 1:N-1
    A[i, i] = -2κ*(1/Δx^2)
end

for i = 1:N-2
    A[i+1, i] = 1κ*(1/Δx^2)
    A[i,i+1] = 1κ*(1/Δx^2)
end

"""
    F(t, x)

Source function as outlined in 3b.
"""
function F(t, x)
    return 2 *(π^2 - 1) * ℯ^(-2t) * sin(π*x) # user specified source function
end

"""
    f(x)

User function for the initial heat distribution
"""
function f(x)
    return sin(π*x) # user specifed initial heat distribution
end

"""
    b(t)

Vector function for user specificed source function.
Passed as argument to euler forward linear.
"""
function b(t)
    return F.(t, x_int)
end

"""
    exact_sol(t, x)

Function for the exact solution.
"""
function exact_sol(t, x)
    return ℯ^(-2t) * sin.(π*x)
end


# Solve y' = Ay + b(t) with y0 = f(x_int)
y0 = f.(x_int)

y = Matrix{Float64}(undef,N+1, M+1)  # my entire solution at all nodes and all time steps.

# fill in initial data:
y[:, 1] = f.(x)

t0 = 0

# Replace below with a call to my_forward_euler_linear()
# forward Euler: y_n+1 = y_n + Δt*(A*y_n + b(t_n))
(t, y) = my_forward_euler_linear(t0, Tf,round(10^-1, digits=5), y0, A, b)

#for n = 1:M
#    # Fill in boundary condition
#    y[1, n+1] = 0
#    y[N+1, n+1] = 0
#    # update remaining rows
#    y[2:N, n+1] = y[2:N, n] + Δt*(A*y[2:N, n] + b(t[n]))
#end
y_actual = zeros(N+1, M+1)
for i = 1:M+1
    y_actual[2:N, i] = exact_sol(t[i], x)[2:N]
end

for i = 1:M+1
    p = plot(x_int, y[:, i], legend = false, ylims=(0,1.5)) # plot the numerical solution
    a = plot!(x_int, y_actual[2:N, i], legend = false, ylims=(0,1.5)) # plot the actual solution
    display(p)
    display(a)

    sleep(0.0002)
    savefig("heatPlot1.png")
end
