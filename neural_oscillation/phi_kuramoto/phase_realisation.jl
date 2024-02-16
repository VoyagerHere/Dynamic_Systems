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

n = [4, 4];
d = [0]
delta = 0.01;
num = length(n);
g_init = 1.01;
alpha = 0.0
tspan = (0, 1000.0)

y0 = zeros(num)
g = [g_init + delta * i for i in  0:(num -  1)];

prob = ODEProblem(eqn!, y0, tspan, (d, alpha, g, n, num))
sol = solve(prob, Tsit5(), reltol=1e-13, abstol=1e-14)

y0 = sol[end]

prob = ODEProblem(eqn!, y0, tspan, (d, alpha, g, n, num))
sol = solve(prob, Tsit5(), reltol=1e-13, abstol=1e-14)

T = sol.t
Y = [mod.(y, 2 * pi) for y in sol.u]
Y = sol.u;
YY = reduce(vcat,transpose.(Y))
YY[:,1]
findall(x -> x < 0, YY[:,1])
plot()

for i in 1:num
  g_cur = g[i]
  n_cur = n[i]
  plot!(T, getindex.(Y, i), label=L"n_%$i = %$n_cur , \gamma_{%$i}=%$g_cur")
end


ylims!(0,  2*pi)
xlabel!(L"t")
ylabel!(L"\varphi")