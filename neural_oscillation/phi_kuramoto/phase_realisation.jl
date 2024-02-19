using DifferentialEquations
using Plots
using LaTeXStrings

Plots.scalefontsizes()
Plots.scalefontsizes(1.5)

function eqn!(du, u, p, t)
    d, alpha, g, n = p
    f = g .- sin.(u ./ n)
    exch = [d * sin(u[2] - u[1] - alpha), d * sin(u[1] - u[2] - alpha)]
    du .= f + exch
end

function eqn!(du, u, p, t)
  d, alpha, g, n, dim_size = p
  f = g .- sin.(u ./ n)
  exch = zeros(dim_size)
  for i in 1:dim_size
    for j in 1:dim_size-1
      exch[i] += d[j] * sin(u[j] - u[i] - alpha)
    end
  end
  du .= f + exch
end

# in 5,5 is process
n = [5, 5];
d = [0.01]
delta = 0.01;
num = length(n);
g_init = 1.01;
alpha = pi/8;
tspan = (1000, 4000)

y0 = zeros(num)
g = [g_init + delta * i for i in  0:(num -  1)];
g = [1.01, 1.02]
prob = ODEProblem(eqn!, y0, tspan, (d, alpha, g, n, num))
sol = solve(prob, Tsit5(), reltol=1e-13, abstol=1e-14)

T = sol.t
# Y = [mod.(y, 2 * pi) for y in sol.u]
n1 = n[1]; g1 = g[1];
plot(T, getindex.(Y, 1), label=L"n_1 = %$n1, \gamma_{1}=%$g1")
for i in 2:num
  n_cur = n[i];
  g_cur = g[i];
  plot!(T, getindex.(Y, i), label=L"n_%$i = %$n_cur , \gamma_{%$i}=%$g_cur")
end

ylims!(0,  2*pi)
xlabel!(L"t")
ylabel!(L"\varphi")