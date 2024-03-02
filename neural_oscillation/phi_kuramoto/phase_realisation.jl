using DifferentialEquations
using Plots
using LaTeXStrings

Plots.scalefontsizes()
Plots.scalefontsizes(1.5)

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

function DELETE_TRANSIENT(Y, tol=0.002)
  len = length(Y);
  i = 10;
  while i < len
      rel_change_1 = abs((Y[i][1] - Y[i-1][1]) / Y[i-1][1])
      rel_change_2 = abs((Y[i][2] - Y[i-1][2]) / Y[i-1][2])
      if ((rel_change_1 < tol) && (rel_change_2 < tol))
          return i
      end
      i += 10;
  end
  return 1
end

n = [2, 2];
d = [0.00]
delta = 0.01;
num = length(n);
g_init = 1.01;
alpha = pi/8;
tspan = (8000, 8500)

y0 = [0, 0]

g = [g_init + delta * i for i in  0:(num -  1)];
g = [1.01, 1.02]
prob = ODEProblem(eqn!, y0, tspan, (d, alpha, g, n, num))
saveat_points = a:0.001:b
sol = solve(prob, Tsit5(), reltol=1e-14, abstol=1e-14, saveat=saveat_points)

T = sol.t
Y = sol.u;

index = DELETE_TRANSIENT(Y)
start = T[index];
y0 = [start, start]

prob = ODEProblem(eqn!, y0, tspan, (d, alpha, g, n, num))
saveat_points = a:0.001:b
sol = solve(prob, Tsit5(), reltol=1e-14, abstol=1e-14, saveat=saveat_points)
Y = sol.u;
T = sol.t;

Y = [mod.(y, 2 * pi) for y in Y]

n1 = n[1]; g1 = g[1];
plot(T, getindex.(Y, 1), label=L"n_1 = %$n1, \gamma_{1}=%$g1")
for i in 2:num
  n_cur = n[i];
  g_cur = g[i];
  plot!(T, getindex.(Y, i), label=L"n_%$i = %$n_cur , \gamma_{%$i}=%$g_cur")
end

# ylims!(0,  2*pi)
xlabel!(L"t")
ylabel!(L"\varphi")