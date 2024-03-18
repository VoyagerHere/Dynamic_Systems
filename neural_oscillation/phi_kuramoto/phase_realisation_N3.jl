using DifferentialEquations
using Plots
using LaTeXStrings

Plots.scalefontsizes()
Plots.scalefontsizes(1.5)

const k_DELETE_TRANSIENT = true; 


function eqn!(dy, y, p, t)
  d, alpha, g, n, dim_size = p
  f = g .- sin.(y ./ n)
  exch = zeros(dim_size)
  d1 = d[1]
  d2 = d[2]
  exch = [d1 * sin(y[2] - y[1] - alpha) + d2 * sin(y[3] - y[1] - alpha); d1 * sin(y[1] - y[2] - alpha) + d2 * sin(y[3] - y[2] - alpha); d1 * sin(y[1] - y[3] - alpha) + d2 * sin(y[2] - y[3] - alpha)];
  dy .= f + exch
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

n = [2, 2, 2];
d1 = 3.0
d = [d1, d1]
num = length(n);
g_init = 1.01;
alpha = 2*pi/3;
tspan = (8000, 8020)

y0 = [0, 0, 0]
const DELTA = 0.005;
G1 = 1.01;
G2 = G1 + DELTA;
G3 = G2 + DELTA

g = [G1, G2, round(G3; digits = 3)]

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
d1 = d[1]
d2 = d[2]
title!(L"\alpha = %$alpha_txt, d_{1} = %$d1, d_{2} = %$d2")
xlabel!(L"t")
ylabel!(L"\varphi")