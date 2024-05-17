using Plots



function eqn_sin(y, t, d, no, F)
  f = g - sin.(y ./ no)
  exch = d * [F[2], F[1]]
  dy_dt = f - exch
  return dy_dt
end


function eqn_cos(y, t, d, no, F)
  f = g - cos.(y ./ no)
  exch = d * [F[2], F[1]]
  dy_dt = f - exch
  return dy_dt
end

eqn = eqn_cos;

function chech_condition(y, sigma)
  y[1] = mod.(y[1], 2 * pi)
  y[2] = mod.(y[2], 2 * pi)


  if ((y[1] > pi/2 - sigma) && (y[1] < pi/2 + sigma))
    F1 =  0;
  else
    F1 = 1;
  end

  if ((y[2] > pi/2 - sigma) && (y[2] < pi/2 + sigma))
    F2 =  0;
  else
    F2 = 1;
  end

  return F1, F2
end

function solver(a, b, sigma, d, y0)
  h = 1/100
  t = a:h:(b-h)
  global F = [0, 0]
  n = length(t)
  y = zeros((n, length(y0)))
  y[1,:] = y0

  for i in 1:n-1
    h = t[i+1] - t[i]
    k1 = eqn(y[i,:], t[i], d, PAR_N, F)
    k2 = eqn(y[i,:] + k1 * h/2, t[i] + h/2, d, PAR_N, F)
    k3 = eqn(y[i,:] + k2 * h/2, t[i] + h/2, d, PAR_N, F)
    k4 = eqn(y[i,:] + k3 * h, t[i] + h, d, PAR_N, F)
    y[i+1,:] = y[i,:] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    F = chech_condition(y[i+1,:], sigma);
  end
  return y, t
end

function draw(T, Y, d, sigma)
  p = plot(T, Y[:, 1], line = (:red), label = "\$\\phi_1\$")
  plot!(T, Y[:, 2], line = (:darkgreen), label = "\$\\phi_2\$")
  ylims!(0, 2 * pi)
  yticks!([0, pi, 2 * pi], ["0", "\\pi", "2\\pi"])
  xlabel!("t")
  ylabel!("\$\\phi\$")
  title!("d = $d, σ = $sigma")
  display(p)
end


d =  0.0005
sigma = 1/2
n1 = 1
n2 = 1
const G1 = 1.001

# DELTA = 0.000005
# DELTA = 0.05
DELTA = 0.5


const G2 = G1+DELTA
N1 = 1 
N2 = 1

const PAR_N = [N1, N2];
g = [G1, G2]
y0 = [pi/2; pi/2]

Y, T = solver(0, 1000, sigma, d, y0)

Y = mod.(Y, 2*pi)
draw(T, Y, d, sigma)