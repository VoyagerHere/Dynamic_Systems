using DifferentialEquations
using Plots
using LaTeXStrings
using JLD
using Statistics
using Dates


Plots.scalefontsizes()
Plots.scalefontsizes(1.5)

const k_ENABLE_ADAPTIVE_GRID = false;
const k_DEBUG_PRINT = false
const k_DRAW_PHASE_REALISATION = false;
const k_IS_SAVE_DATA = true;
const k_DELETE_TRANSIENT = false;
const k_DELETE_UNSTABLE = false;

const DATA_TAKE_ERROR = 0.05;

global a = 8000;
global b = 9000;

# For ADAPTIVE_GRID
const init_b = 9000;
const b_step = 3000;

const ADAPTIVE_SET_ERROR = 10;
const SPIKE_ERROR =  0


name = "pi_8__2_2"
const ALPHA = pi / 8
N1 = 2
N2 = 2

const NUM = 2;
global PAR_N = [N1, N2];
const D_MAX =  0.04
const D_ACCURACY =  0.00001
const G_NUM = 500
const SYNC_ERROR =  0.05
const GStart =  1.01
const DELTA = 0.005;
D_LIST = 0:D_ACCURACY:D_MAX
D_NUM = length(D_LIST)

G1 = 1.01;
G2 = G1 + DELTA;

const NUM_OF_COMPUTE_RES = 3;
DATA = [zeros(NUM_OF_COMPUTE_RES) for _ in 1:D_NUM]
W = [zeros(NUM*2) for _ in 1:D_NUM]

const ALPHA_TEXT = L"π/8"

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

function FREQ_SYNC(DATA, G1, G2, PAR_N, NUM, D_LIST, SPIKE_ERROR, ALPHA)
  num_of_iterations = length(D_LIST)
  k_ENABLE_ADAPTIVE_GRID && ADAPTIVE_GRID(0, true);
  for m in eachindex(D_LIST)
    global D = D_LIST[m]
    
    tspan = (a, b)
    
    p = (D, ALPHA, [G1, G2], PAR_N, NUM);
    y0 = [0; 0]

    prob = ODEProblem(eqn!, y0, tspan, p)
    saveat_points = a:0.001:b
    sol = solve(prob, Tsit5(), reltol=1e-14, abstol=1e-14, saveat=saveat_points)
    global Y = sol.u;
    global T = sol.t;

    if (k_DELETE_TRANSIENT)
      index = DELETE_TRANSIENT(Y)
      start = T[index];
      y0 = [start, start]

      prob = ODEProblem(eqn!, y0, tspan, p)
      sol = solve(prob, Tsit5(), reltol=1e-14, abstol=1e-14)
      Y = sol.u;
      T = sol.t;
    end

    w_s, w_b, err = SYNC_PAIR(T, Y, PAR_N, SPIKE_ERROR)
    ratio_s = w_s[1] / w_s[2];
    ratio_b = w_b[1] / w_b[2];

    if (err != 0)
      DATA[m] = [D, -err, -err]
      continue;
    end

    DATA[m] = [D, ratio_s, ratio_b]
    W[m] = vcat(w_s, w_b)
    println("Iteration $m of $num_of_iterations")
    return;
  end
end

function SYNC_PAIR(T, Y, PAR_N, error)
  Y = reduce(vcat,transpose.(Y))
  SPIKES1, err1 = FIND_SPIKES(Y[:,1], PAR_N[1])
  SPIKES2, err2 = FIND_SPIKES(Y[:,2], PAR_N[2])

  if (k_DELETE_UNSTABLE)
    unstbl_1 = DELETE_UNSTBL(SPIKES1, err1, PAR_N[1], error)
    unstbl_2 = DELETE_UNSTBL(SPIKES1, err2, PAR_N[2], error)

    untsbl_ind = max(unstbl_1, unstbl_2)

    SPIKES1 = SPIKES1[untsbl_ind:end];
    SPIKES2 = SPIKES2[untsbl_ind:end];
  end

  k_DEBUG_PRINT && println("Spikes in 1: ", length(SPIKES1))
  k_DEBUG_PRINT && println("Spikes in 2: ", length(SPIKES2))

  err = handle_errors(Bool(err1), Bool(err2));
  if err != 0
      return ([0, 0], [0, 0], err)
  end

  k_ENABLE_ADAPTIVE_GRID && (ADAPTIVE_GRID(minimum([length(SPIKES1), length(SPIKES1)]), false))

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

function ADAPTIVE_GRID(num_of_bursts, reset)
  global b;
  if (reset == true)
    b = init_b;
  else   
    if (num_of_bursts < ADAPTIVE_SET_ERROR)
      b += b_step;
    end
  end
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

function DELETE_UNSTBL(SPIKES, err, n, unstable_err)
  UNSTBL = 1
  if err > 0
      return UNSTBL
  end
  B = zeros(Int64, div(length(SPIKES), n))
  for m in n:n:length(SPIKES)
    B[div(m, n)] = SPIKES[m]
  end
  if length(B) < unstable_err + 2
      return UNSTBL
  end
  element = B[unstable_err]
  UNSTBL = findall(SPIKES .== element)
  if isempty(UNSTBL)
      UNSTBL = 1
  end
  return UNSTBL[1]
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
  plot(T, getindex.(Y, 1), label=L"n_1 = %$n1, \gamma_{1}=%$G1")
  plot!(T, getindex.(Y, 2), label=L"n_2 = %$n2 , \gamma_{2}=%$G2")
  title!(L"d = %$D")
  ylims!(0,  2*pi)
  xlabel!(L"t")
  ylabel!(L"\varphi")
end

FREQ_SYNC(DATA, G1, G2, PAR_N, NUM, D_LIST, SPIKE_ERROR, ALPHA);


if k_IS_SAVE_DATA 
  times = Dates.format(now(),"__yyyymmdd_HHMM");
  filename ="$name$times.jld2"
  @save filename DATA W
end