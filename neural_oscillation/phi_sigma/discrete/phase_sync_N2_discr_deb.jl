using DifferentialEquations
using Plots
using LaTeXStrings
using JLD
using Statistics
using Dates


const k_ENABLE_ADAPTIVE_GRID = false;
const k_DEBUG_PRINT = false
const k_DRAW_PHASE_REALISATION = false;
const k_IS_SAVE_DATA = true;
const k_DELETE_TRANSIENT = true; 
const k_PRINT_ITERATION = false;

const DATA_TAKE_ERROR = 0.25;

# For ADAPTIVE_GRID
const b_step = 2000;
const ADAPTIVE_SET_ERROR = 10;

# For DELETE_UNSTABLE
const SPIKE_ERROR = 0

N1 = 1 
N2 = 1

const G1 = 1.001

# DELTA = 0.005
DELTA = 0.01
# DELTA = 0.05

const G2 = G1+DELTA


func_txt = "sin"

name = "ph_discr_$func_txt$N1$N2$DELTA"

g = [G1, G2]
const D_ACCURACY =  0.001
const sigma_ACCURACY =  0.001
const D_MAX =  0.1
const sigma_MAX = 1

const NUM = 2;
const PAR_N = [N1, N2];
const SYNC_ERROR =  0.25;
sigma_LIST = 0.812:sigma_ACCURACY:sigma_MAX
D_LIST = 0.003:D_ACCURACY:D_MAX
D_NUM = length(D_LIST)
sigma_NUM = length(sigma_LIST)

const NUM_OF_COMPUTE_RES = 3;
DATA = [zeros(NUM_OF_COMPUTE_RES) for _ in 1:(D_NUM*sigma_NUM)]
SYNC = [zeros(2) for _ in 1:(D_NUM*sigma_NUM)]


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

eqn = eqn_sin;

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

function PHASE_SYNC(DATA, SYNC, PAR_N, sigma_LIST, D_LIST, SPIKE_ERROR)
  for m in eachindex(D_LIST)
      global D = D_LIST[m]
      a = 8000;
      b = 10000;
      for k in eachindex(sigma_LIST)
        sigma = sigma_LIST[k];
        y0 = [pi/2; pi/2]
        global Y, T = solver(a, b, sigma, D, y0)
        if (k_DELETE_TRANSIENT)
          y0 = Y[end, :]
          global Y, T = solver(a, b, sigma, D, y0)
        end 

        global DIFF_SP, DIFF_BS, ratio, err, b_new, len = SYNC_PAIR(T, Y, PAR_N, SPIKE_ERROR, b)

        b = copy(b_new);
        sync = [0, 0]
        if (err == 0)
          sync[1] = IS_SYNC(DIFF_BS, SYNC_ERROR, ratio);
          sync[2] = IS_SYNC(DIFF_SP, SYNC_ERROR, ratio);
          if (sum(sync) == 0)
            ratio = 0
          end
        end

        @inbounds DATA[m + (k-1)*D_NUM] = [D, ratio, sigma]
        @inbounds SYNC[m + (k-1)*D_NUM] = sync;
        return;
    end
    k_PRINT_ITERATION && println("Iteration $k of $num_of_iterations")
    end
end

function SYNC_PAIR(T, Y, PAR_N, error, b)
  SPIKES1, err1 = FIND_SPIKES(Y[:,1], PAR_N[1])
  SPIKES2, err2 = FIND_SPIKES(Y[:,2], PAR_N[2])


  len = [length(SPIKES1), length(SPIKES2)];

  err = handle_errors(Bool(err1), Bool(err2));
  if err != 0
      DIFF_SP = 0
      DIFF_BS = 0
      ratio = -err
      return (DIFF_SP, DIFF_BS, ratio, err, b, len)
  end

  BURSTS1 = FIND_BURST(SPIKES1, PAR_N[1])
  BURSTS2 = FIND_BURST(SPIKES2, PAR_N[2])
  
  k_DEBUG_PRINT && println("Bursts in 1: ", length(BURSTS1))
  k_DEBUG_PRINT && println("Bursts in 2: ", length(BURSTS2))

  k_ENABLE_ADAPTIVE_GRID && (b = ADAPTIVE_GRID(minimum([length(BURSTS1), length(BURSTS2)]), b))
  
  ratio = FIND_RATIO(BURSTS1, BURSTS2)

  DIFF_SP = FIND_DIFF(SPIKES1, SPIKES2, T)
  DIFF_BS = FIND_DIFF(BURSTS1, BURSTS2, T)
  return (DIFF_SP, DIFF_BS, ratio, err, b, len)
end

function ADAPTIVE_GRID(num_of_bursts, b)
  if (num_of_bursts < (ADAPTIVE_SET_ERROR + SPIKE_ERROR))
     b = b + b_step;
     return b;
  else
    return b;
  end
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

function FIND_BURST(SPIKES, n)
  B = zeros(Int64, div(length(SPIKES), n))
  for m in n:n:length(SPIKES)
      B[div(m, n)] = SPIKES[m]
  end
  return B
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
  if length(B) < 2*unstable_err 
    # err = 1;
    return UNSTBL
  end
  element = B[unstable_err]
  UNSTBL = findall(SPIKES .== element)
  if isempty(UNSTBL)
      UNSTBL = 1
  end
  return UNSTBL[1]
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

function FIND_RATIO(A, B)
    ratio = zeros(Int64, length(A) - 2)
    for i in 2:length(A) - 1
        ratio[i - 1] = sum(B .< A[i + 1]) - sum(B .< A[i])
    end
    return Int64.(floor(mean(ratio)))
end

# Find near spike by left for spike from A1
# and compute diff as difference between them
# A2 have bigger gamma value
function FIND_DIFF(A1, A2, T)
    NEAR = []
    if (minimum.(A1) > minimum.(A2)) 
      A1, A2 = A2, A1
    end

    for el in A2
        _, index = findmax(A1[A1 .<= el])
        if !isempty(index)
            push!(NEAR, A1[index])
        end
    end
    A2 = T[A2];
    NEAR = T[NEAR];

    len = length(NEAR)
    DIFF = zeros(len)

    # A2 = A2[end - len + 1:end]
    DIFF = A2 .- NEAR
    return DIFF
end

function IS_SYNC(DIFF, SYNC_ERROR, ratio)
  if ((ratio == 1) &&  DIFF_SPIKES(DIFF, SYNC_ERROR))
    return 2;
  end
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


@time  PHASE_SYNC(DATA, SYNC, PAR_N, sigma_LIST, D_LIST, SPIKE_ERROR)

if k_IS_SAVE_DATA 
  times = Dates.format(now(),"__yyyymmdd_HHMM");
  filename ="$name$times.jld2"
  @save filename DATA SYNC
end