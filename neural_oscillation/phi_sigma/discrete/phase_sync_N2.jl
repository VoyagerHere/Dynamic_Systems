using DifferentialEquations
using Plots
using LaTeXStrings
using JLD
using Statistics
using Dates


Plots.scalefontsizes()
Plots.scalefontsizes(1.5)

const k_ENABLE_ADAPTIVE_GRID = true;
const k_DEBUG_PRINT = false
const k_DRAW_PHASE_REALISATION = false;
const k_IS_SAVE_DATA = true;
const k_DELETE_TRANSIENT = false; 
const k_DELETE_UNSTABLE = true;

const DATA_TAKE_ERROR = 0.05;

global a = 8000;
global b = 10000;

# For ADAPTIVE_GRID #обрезает берсты
const init_b = 10000;
const b_step = 3000;
const ADAPTIVE_SET_ERROR = 10;

# For DELETE_UNSTABLE
const SPIKE_ERROR =  10

name = "untitled"
N1 = 2 
N2 = 2
const NUM = 2;
const K = -500
global PAR_N = [N1, N2];
const D_MAX =  1.5
const D_ACCURACY =  0.005
DELTA = 0.001
const gamma1 = 1.01
const gamma2 = 1.01+DELTA
const SYNC_ERROR =  0.05
const sigma_MAX = pi
const sigma_ACCURACY =  0.005
sigma_LIST = 0:sigma_ACCURACY:sigma_MAX
D_LIST = 0:D_ACCURACY:D_MAX
D_NUM = length(D_LIST)
sigma_NUM = length(sigma_LIST)

const NUM_OF_COMPUTE_RES = 3;
DATA = [zeros(NUM_OF_COMPUTE_RES) for _ in 1:(D_NUM*sigma_NUM)]
SYNC = [zeros(2) for _ in 1:(D_NUM*sigma_NUM)]


function eqn!(du, u, p, t)
  d, sigma, g, n, dim_size = p
  F1, F2 = check_condition(u, sigma)
  f = g .- sin.(u ./ n)
  exch = d*[F2, F1]
  du .= f - exch
end


function check_condition(y, sigma)
  F1 = y[1] > pi/2 - sigma && y[1] < pi/2 + sigma ?  0 :  1
  F2 = y[2] > pi/2 - sigma && y[2] < pi/2 + sigma ?  0 :  1
  return F1, F2
end

function PHASE_SYNC(DATA, SYNC, PAR_N, NUM, sigma_LIST, D_LIST, SPIKE_ERROR)
    num_of_iterations = length(sigma_LIST)
    for k in eachindex(sigma_LIST)
      sigma = sigma_LIST[k];
      k_ENABLE_ADAPTIVE_GRID && ADAPTIVE_GRID(0, true);
      for m in eachindex(D_LIST)
        d = D_LIST[m]
        
        tspan = (a, b)
        
        p = (d, sigma, [gamma1, gamma2],  PAR_N, NUM);
        y0 = [pi/2; pi/2]

        prob = ODEProblem(eqn!, y0, tspan, p)
        sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
        Y = sol.u;
        T = sol.t;

        if (k_DELETE_TRANSIENT)
          index = DELETE_TRANSIENT(Y)
          start = T[index];
          y0 = [start, start]

          prob = ODEProblem(eqn!, y0, tspan, p)
          sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
          Y = sol.u;
          T = sol.t;
        end 

        DIFF_SP, DIFF_BS, ratio, err = SYNC_PAIR(T, Y, PAR_N, SPIKE_ERROR)

        sync = [0, 0]
        if (err == 0)
          sync[1] = IS_SYNC(DIFF_BS, SYNC_ERROR);
          sync[2] = IS_SYNC(DIFF_SP, SYNC_ERROR);
          if (sum(sync) == 0)
            ratio = 0
          end
        end
        
        
        DATA[m + (k-1)*D_NUM] = [d, ratio, sigma]
        SYNC[m + (k-1)*D_NUM] = sync;
    end
      println("Iteration $k of $num_of_iterations")
    end
end

function SYNC_PAIR(T, Y, PAR_N, error)
  Y = reduce(vcat,transpose.(Y))
  SPIKES1, err1 = FIND_SPIKES(Y[:,1], PAR_N[1])
  SPIKES2, err2 = FIND_SPIKES(Y[:,2], PAR_N[2])

  if (k_DELETE_UNSTABLE)
    unstbl_1 = DELETE_UNSTBL(SPIKES1, err1, PAR_N[1], error)
    unstbl_2 = DELETE_UNSTBL(SPIKES2, err2, PAR_N[2], error)

    SPIKES1 = SPIKES1[unstbl_1:end];
    SPIKES2 = SPIKES2[unstbl_2:end];
  end

  err = handle_errors(Bool(err1), Bool(err2));
  if err != 0
      DIFF_SP = 0
      DIFF_BS = 0
      ratio = -err
      return (DIFF_SP, DIFF_BS, ratio, err)
  end

  BURSTS1 = FIND_BURST(SPIKES1, PAR_N[1])
  BURSTS2 = FIND_BURST(SPIKES2, PAR_N[2])
  k_DEBUG_PRINT && println("Bursts in 1: ", length(BURSTS1))
  k_DEBUG_PRINT && println("Bursts in 2: ", length(BURSTS2))

  k_ENABLE_ADAPTIVE_GRID && (ADAPTIVE_GRID(minimum([length(BURSTS1), length(BURSTS2)]), false))
  
  ratio = FIND_RATIO(BURSTS1, BURSTS2)

  DIFF_SP = FIND_DIFF(SPIKES1, SPIKES2, T)
  DIFF_BS = FIND_DIFF(BURSTS1, BURSTS2, T)
  return (DIFF_SP, DIFF_BS, ratio, err)
end

function ADAPTIVE_GRID(num_of_bursts, reset)
  global b;
  if (reset == true)
    b = init_b;
  else   
    if (num_of_bursts < (ADAPTIVE_SET_ERROR + SPIKE_ERROR))
      b += b_step;
    end
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

function DRAW(T, Y, gamma1, gamma2, D, PAR_N)
    Y = [mod.(y, 2 * pi) for y in Y]
    n1 = PAR_N[1];
    n2 = PAR_N[2];
    plot(T, getindex.(Y, 1), label=L"n_1 = %$n1, \gamma_{1}=%$gamma1")
    plot!(T, getindex.(Y, 2), label=L"n_2 = %$n2 , \gamma_{2}=%$gamma2")
    title!(L"d = %$D")
    ylims!(0,  2*pi)
    xlabel!(L"t")
    ylabel!(L"\varphi")
end

PHASE_SYNC(DATA, SYNC, PAR_N, NUM, sigma_LIST, D_LIST, SPIKE_ERROR)

if k_IS_SAVE_DATA 
  times = Dates.format(now(),"__yyyymmdd_HHMM");
  filename ="$name$times.jld2"
  @save filename DATA SYNC
end