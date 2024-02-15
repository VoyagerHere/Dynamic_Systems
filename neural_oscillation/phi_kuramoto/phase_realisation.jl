# solver = CVODE_BDF()

# # Set solver options
# solver_options = Dict{Symbol, Any}(
#     :reltol =>  1e-8,
#     :abstol =>  1e-8,
#     :maxstep =>  0.1
# )

# # Solve the problem with the chosen solver and options
# sol = solve(prob, solver, solver_options)
import Pkg; Pkg.add("DifferentialEquations")
using DifferentialEquations
using Plots


function eqn(dy, y, p, t)
  g, d, alpha, n = p
  f = g .- sin.(y ./ n)
  exch = [d * sin(y[2] - y[1] - alpha), d * sin(y[1] - y[2] - alpha)]
  dy .= f .+ exch
end

# Initial conditions
y0 = [0.0,  0.0]

# Parameters
d =  0.08
delta =  0.01
alpha =  2*pi/3
n = [3,  3]

p = (g1, delta, alpha, n)
tspan = (2000,  4000)
prob = ODEProblem(eqn, y0, tspan, p)

sol = solve(prob, Tsit5())

Y = mod.(sol.u,  2*pi)

plot(sol.t, Y[:,1], label="\$\\phi_1, n_1 =  3, \\gamma_1 =  1.01\$", lw=2)
plot!(sol.t, Y[:,2], label="\$\\phi_2, n_2 =  3, \\gamma_2 =  1.02\$", lw=2)
xlabel!("t")
ylabel!("phi")
title!("Delta = $(delta), d = $(d), alpha = $alpha")
ylims!(0,  2*pi)
yticks!([0, pi,  2*pi], ["0", "\$\\pi\$", "2\$\\pi\$"])
legend!(loc=:upperleft, frame=:box)
