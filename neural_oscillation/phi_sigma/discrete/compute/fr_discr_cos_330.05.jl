using DifferentialEquations
using Plots
using LaTeXStrings
using JLD
using Statistics
using Dates


const k_DEBUG_PRINT = false
const k_DRAW_PHASE_REALISATION = false;
const k_IS_SAVE_DATA = true;
const k_DELETE_TRANSIENT = true; 
const k_PRINT_ITERATION = false;

const DATA_TAKE_ERROR = 0.25;


N1 = 3
N2 = 3

# DELTA = 0.005
# DELTA = 0.01
DELTA = 0.05


SIGMA_FIXED = 1/2;

const D_MAX =  0.03

func_txt = "cos"

name = "fr_dicr_$func_txt$N1$N2$DELTA"


const G1 = 1.001
const G2 = G1 + DELTA
g = [G1, G2]
const D_ACCURACY =  0.0001


const NUM = 2;
const PAR_N = [N1, N2];
const SYNC_ERROR = 0.25;
D_LIST = 0:D_ACCURACY:D_MAX
D_NUM = length(D_LIST)


DATA = [zeros(4) for _ in 1:(D_NUM)]
W = [zeros(2) for _ in 1:(D_NUM)]

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

function chech_condition(y, sigma, PAR_N)
  y[1] = mod.(y[1], 2 * pi)
  y[2] = mod.(y[2], 2 * pi)

  y = y ./ PAR_N

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
    F = chech_condition(y[i+1,:], sigma, PAR_N);
  end
  return y, t
end

function FREQ_SYNC(DATA, PAR_N, D_LIST)

    sigma = SIGMA_FIXED;
    a = 8000;
    b = 12000;

    for m in eachindex(D_LIST)    
      D = D_LIST[m]
      println(D)
      y0 = [pi/2; pi/2]
      Y, T = solver(a, b, sigma, D, y0)

      if (k_DELETE_TRANSIENT)
        y0 = Y[end, :]
        Y, T = solver(a, b, sigma, D, y0)
      end 

      w_s, w_b, err = SYNC_PAIR(T, Y, PAR_N)

      ratio_s = w_s[1] / w_s[2];
      ratio_b = w_b[1] / w_b[2];

      if (err != 0)
        DATA[m] = [D,  -err, -err]
        continue;
      end
        DATA[m] = [D, ratio_s, ratio_b]
        W[m] = vcat(w_s, w_b)
      end
    
end

function SYNC_PAIR(T, Y, PAR_N)
  global SPIKES1, err1 = FIND_SPIKES(Y[:,1], PAR_N[1])
  global SPIKES2, err2 = FIND_SPIKES(Y[:,2], PAR_N[2])

  k_DEBUG_PRINT && println("Spikes in 1: ", length(SPIKES1))
  k_DEBUG_PRINT && println("Spikes in 2: ", length(SPIKES2))

  err = handle_errors(Bool(err1), Bool(err2));
  if err != 0
      return ([0, 0], [0, 0], err)
  end

  Times_1 = FIND_TIMES(SPIKES1, T, PAR_N);
  Times_2 = FIND_TIMES(SPIKES2, T, PAR_N);

  Times_SP_1 = DEL_MIDDLE_BURST_INTRVL(PAR_N[1], Times_1)
  Times_SP_2 = DEL_MIDDLE_BURST_INTRVL(PAR_N[2], Times_2)

  avg(x) = (ones(length(x)) / length(x))'*x

  wb_1 = 2*pi./sum(Times_1);
  wb_2 = 2*pi./sum(Times_2);

  ws_1 = 2*pi./avg(Times_SP_1);
  ws_2 = 2*pi./avg(Times_SP_2);

  return ([ws_1, ws_2], [wb_1, wb_2], err)
end

avg(x) = (ones(length(x)) / length(x))'*x

function FIND_SPIKES(Y, n)
  Y = mod.(Y, 2*pi)
  SPIKES = findall(x -> abs.(x - 2*pi) < DATA_TAKE_ERROR, Y )
  FIND_NEAR_POINTS(SPIKES)
  if (length(SPIKES)/n < 3)
    err = 1;
  else
    err = 0;
  end
  return (SPIKES, err)
end


function FIND_TIMES(SPIKES, T, NUM)
  len = maximum(NUM);
  TIMES = zeros(len);
  for k in 2:(len+1)
    val1 = SPIKES[k-1]
    val2 = SPIKES[k]
    TIMES[k-1] = T[val2] - T[val1];
  end
  return TIMES;
end

function DEL_MIDDLE_BURST_INTRVL(n, T)
  if (n > 1)
    T = sort(T, rev = true)
    T = T[2:end]
  end
  return T;
end

function handle_errors(err1::Bool, err2::Bool)
  if err1 && err2
      return  3
  elseif err2
      return  2
  elseif err1
      return  1
  else
      return  0
  end
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

function FIND_NEAR_POINTS(POINTS)
  i = 1;
    while i < length(POINTS)
        if POINTS[i+1] - POINTS[i] == 1
            deleteat!(POINTS, i)
        else
            i +=  1
        end
    end
end

function IS_SYNC(DIFF, SYNC_ERROR)
  if (DIFF_SPIKES(DIFF, SYNC_ERROR))
    return 1;
  end
  return 0;
end

function DIFF_SPIKES(DIFF, SYNC_ERROR)
  mn = mean(abs.(DIFF))
  return all(abs.(abs.(DIFF) .- mn) .< SYNC_ERROR)
end

function DRAW(T, Y, G1, G2, D, PAR_N)
  Y = [mod.(y, 2 * pi) for y in Y]
  n1 = PAR_N[1];
  n2 = PAR_N[2];
  plot(T, Y[:,2], label=L"n_1 = %$n1, \gamma_{1}=%$G1", linecolor=:red)
  plot!(T, Y[:,1], label=L"n_2 = %$n2 , \gamma_{2}=%$G2", linecolor=:darkgreen)
  title!(L"d = %$D")
  ylims!(0,  2*pi)
  xlabel!(L"t")
  ylabel!(L"\varphi")
end

FREQ_SYNC(DATA, PAR_N, D_LIST);


if k_IS_SAVE_DATA 
  times = Dates.format(now(),"__yyyymmdd_HHMM");
  filename ="$name$times.jld2"
  @save filename DATA W
end