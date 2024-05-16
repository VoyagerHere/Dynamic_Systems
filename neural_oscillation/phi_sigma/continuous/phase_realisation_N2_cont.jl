using DifferentialEquations
using Plots
using LaTeXStrings

Plots.scalefontsizes()
Plots.scalefontsizes(1.5)

const k_DELETE_TRANSIENT = false; 

K = -500


function eqn!(du, u, p, t)#u - это тета
  d, sigma, g, n = p
  f = g .- sin.(u ./ n)
  exch = 1 ./ (1 .+ ℯ.^(K*(cos(sigma).-sin.(u))))
  du .= f - d*exch
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

n = [3, 3];
d = 0.01
sigma = 1/2
tspan = (8000, 9000)

y0 = [0, 0]

# const DELTA =  0.05
# const DELTA = 0.000005
DELTA = 0.000005

G1 = 1.01;
G2 = G1 + DELTA;

g = [G1, G2]

p = (d, sigma, g, n)
prob = ODEProblem(eqn!, y0, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-14, abstol=1e-14)

T = sol.t
Y = sol.u;


if (k_DELETE_TRANSIENT)
  y0 = Y[end];

  prob = ODEProblem(eqn!, y0, tspan, (d, sigma, g, n))
  sol = solve(prob, Tsit5(), reltol=1e-14, abstol=1e-14)
  Y = sol.u;
  T = sol.t;
end

Y = [mod.(y, 2 * pi) for y in Y]

n1 = n[1]; g1 = g[1];
n2 = n[2]; g2 = g[2];

plot(T, getindex.(Y, 1), label=L"n_1 = %$n1, \gamma_{1}=%$g1")
plot!(T, getindex.(Y, 2), label=L"n_2 = %$n2, \gamma_{2}=%$g2")


# ylims!(0,  2*pi)
title!(L"\sigma = %$sigma, d = %$d")
xlabel!(L"t")
ylabel!(L"\varphi")
