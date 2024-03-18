using DifferentialEquations
using Plots
using LaTeXStrings

Plots.scalefontsizes()
Plots.scalefontsizes(1.5)

const k_DELETE_TRANSIENT = true; 


function eqn!(du, u, p, t)
  d, alpha, g, n = p
  f = g .- sin.(u ./ n)
  exch = [d * sin(u[2] - u[1] - alpha), d * sin(u[1] - u[2] - alpha)]
  du .= f .+ exch
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
d = 7.5
num = length(n);
alpha = pi/8;
alpha_txt = "π/8"
tspan = (8030, 8130)

y0 = [0, 0]

const DELTA = 0.005;
G1 = 1.01;
G2 = G1 + DELTA;

g = [G1, G2]
prob = ODEProblem(eqn!, y0, tspan, (d, alpha, g, n, num))
sol = solve(prob, Tsit5(), reltol=1e-14, abstol=1e-14)

T = sol.t
Y = sol.u;


if (k_DELETE_TRANSIENT)
  y0 = Y[end];

  prob = ODEProblem(eqn!, y0, tspan, (d, alpha, g, n, num))
  sol = solve(prob, Tsit5(), reltol=1e-14, abstol=1e-14)
  Y = sol.u;
  T = sol.t;
end

Y = [mod.(y, 2 * pi) for y in Y]

n1 = n[1]; g1 = g[1];
plot(T, getindex.(Y, 1), label=L"n_1 = %$n1, \gamma_{1}=%$g1")
for i in 2:num
  n_cur = n[i];
  g_cur = g[i];
  plot!(T, getindex.(Y, i), label=L"n_%$i = %$n_cur , \gamma_{%$i}=%$g_cur")
end

# ylims!(0,  2*pi)
title!(L"\alpha = %$alpha_txt, d = %$d")
xlabel!(L"t")
ylabel!(L"\varphi")
