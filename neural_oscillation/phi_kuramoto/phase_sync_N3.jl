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
const k_IS_SAVE_DATA = false;
const k_DELETE_TRANSIENT = false;
const k_DELETE_UNSTABLE = false;

const DATA_TAKE_ERROR = 0.05;


global a = 1000;
global b = 2000;

# For ADAPTIVE_GRID
const init_b = 2000;
const b_step = 1000;
const ADAPTIVE_SET_ERROR = 10;

const SPIKE_ERROR =  10


N1 = 1;
N2 = 1;
N3 = 1;
const NUM = 3;
global PAR_N = [N1, N2, N3];
const D_MAX =  0.07
const D_ACCURACY =  0.0001
const G_NUM = 500
const SYNC_ERROR =  0.05
const GStart =  1.01
const DELTA =  0.025
G_LIST = range(GStart, stop=GStart + DELTA, length=G_NUM)
D_LIST = 0:D_ACCURACY:D_MAX
D_NUM = length(D_LIST)

const NUM_OF_COMPUTE_RES = 4;
DATA = [zeros(NUM_OF_COMPUTE_RES) for _ in 1:(D_NUM*G_NUM)]
SYNC = [zeros(5) for _ in 1:(D_NUM*G_NUM)]
DEATH = [zeros(NUM) for _ in 1:(D_NUM*G_NUM)]


const ALPHA_TEXT = L"π/8"
const ALPHA = pi /  8

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

function PHASE_SYNC(DATA, SYNC, GStart, PAR_N, NUM, G_LIST, D_LIST, SPIKE_ERROR, ALPHA)
    num_of_iterations = length(G_LIST)
    G1 = GStart;
    for k in eachindex(G_LIST)
      G2 = G_LIST[k];
      G3 = G2 + DELTA;
      k_ENABLE_ADAPTIVE_GRID && ADAPTIVE_GRID(0, true);
      for m in eachindex(D_LIST)
        d1 = D_LIST[m]
        d2 = d1;
        
        tspan = (a, b)
        
        p = ([d1, d2], ALPHA, [G1, G2, G3], PAR_N, NUM);
        y0 = [0; 0; 0]

        prob = ODEProblem(eqn!, y0, tspan, p)
        sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
        Y = sol.u;
        T = sol.t;

        if (k_DELETE_TRANSIENT)
          index = DELETE_TRANSIENT(Y)
          start = T[index];
          y0 = [start, start, start]

          prob = ODEProblem(eqn!, y0, tspan, p)
          sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
          Y = sol.u;
          T = sol.t;
        end 

        DIFF_SP_12, DIFF_BS_12, ratio_12, err_12 = SYNC_PAIR(T, Y, PAR_N, SPIKE_ERROR, 1, 2)
        DIFF_SP_23, DIFF_BS_23, ratio_23, err_23 = SYNC_PAIR(T, Y, PAR_N, SPIKE_ERROR, 2, 3)

        err = vcat(err_12, err_23);
        deleteat!(err, 2); # remove double record of second neuron

        sync = zeros(NUM * 2);
        if (sum(err_12) == 0)
          sync[1] = IS_SYNC(DIFF_BS_12, SYNC_ERROR);
          sync[2] = IS_SYNC(DIFF_SP_12, SYNC_ERROR);
        end
        if (sum(err_23) == 0)
          sync[3] = IS_SYNC(DIFF_BS_23, SYNC_ERROR);
          sync[4] = IS_SYNC(DIFF_SP_23, SYNC_ERROR);
        end
        if (sum(err) == 0)
          sync[5] = IS_SYNC([DIFF_BS_12; DIFF_BS_23], SYNC_ERROR);
          sync[6] = IS_SYNC([DIFF_SP_12; DIFF_SP_23], SYNC_ERROR);
        end

        delta = G2 - G1
        DATA[m + (k-1)*D_NUM] = [d1, d2, ratio_12, ratio_23, delta]
        SYNC[m + (k-1)*D_NUM] = sync;
        DEATH[m + (k-1)*D_NUM] = err;
    end
      println("Iteration $k of $num_of_iterations")
    end
end

function SYNC_PAIR(T, Y, PAR_N, error, ind1, ind2)
  Y = reduce(vcat,transpose.(Y))
  SPIKES1, err1 = FIND_SPIKES(Y[:,ind1], PAR_N[ind1])
  SPIKES2, err2 = FIND_SPIKES(Y[:,ind2], PAR_N[ind2])

  if (k_DELETE_UNSTABLE)
    unstbl_1 = DELETE_UNSTBL(BURSTS1, err1, PAR_N[1], error)
    unstbl_2 = DELETE_UNSTBL(BURSTS2, err2, PAR_N[2], error)

    SPIKES1 = SPIKES1[unstbl_1:end];
    SPIKES2 = SPIKES2[unstbl_2:end];
  end

  k_DEBUG_PRINT && println("Spikes in $ind1: ", length(SPIKES1))
  k_DEBUG_PRINT && println("Spikes in $ind2: ", length(SPIKES2))

  global err = [err1, err2];

  if sum(err) > 0 
      DIFF_SP = 0
      DIFF_BS = 0
      ratio = -1
      return (DIFF_SP, DIFF_BS, ratio, err)
  end

  global BURSTS1 = FIND_BURST(SPIKES1, PAR_N[ind1])
  global BURSTS2 = FIND_BURST(SPIKES2, PAR_N[ind2])

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
    if (num_of_bursts < ADAPTIVE_SET_ERROR)
      b += b_step;
    end
  end
end

function FIND_SPIKES(Y, n)
  Y = mod.(Y, 2*pi)
  global SPIKES = findall(x -> abs.(x - 2*pi) < DATA_TAKE_ERROR, Y )
  FIND_NEAR_POINTS(SPIKES)
  if (length(SPIKES)/n < 2)
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
  if length(B) < unstable_err + 2
      return UNSTBL
  end
  element = B[unstable_err]
  global UNSTBL = findall(SPIKES .== element)
  if isempty(UNSTBL)
      UNSTBL = 1
  end
  return UNSTBL[1]
end

function DELETE_TRANSIENT(Y, tol=0.002)
  len = length(Y);
  i = 10;
  while i < len
      global rel_change_1 = abs((Y[i][1] - Y[i-1][1]) / Y[i-1][1])
      global rel_change_2 = abs((Y[i][2] - Y[i-1][2]) / Y[i-1][2])
      global rel_change_3 = abs((Y[i][3] - Y[i-1][3]) / Y[i-1][3])
      if ((rel_change_1 < tol) && (rel_change_2 < tol) && (rel_change_3 < tol))
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
    ratio = zeros(length(A) - 2)
    for i in 2:length(A) - 1
        ratio[i - 1] = sum(B .< A[i + 1]) - sum(B .< A[i])
    end
    return mean(ratio)
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

function DRAW(T, Y, G1, G2, G3, D, PAR_N)
  Y = [mod.(y, 2 * pi) for y in Y]
  n1 = PAR_N[1];
  n2 = PAR_N[2];
  n3 = PAR_N[3];
  plot(T, getindex.(Y, 1), label=L"n_1 = %$n1, \gamma_{1}=%$G1")
  plot!(T, getindex.(Y, 2), label=L"n_2 = %$n2 , \gamma_{2}=%$G2")
  plot!(T, getindex.(Y, 3), label=L"n_3 = %$n3 , \gamma_{2}=%$G3")
  title!(L"d = %$D")
  ylims!(0,  2*pi)
  xlabel!(L"t")
  ylabel!(L"\varphi")
end

PHASE_SYNC(DATA, SYNC, GStart, PAR_N, NUM, G_LIST, D_LIST, SPIKE_ERROR, ALPHA);

if k_IS_SAVE_DATA 
  time = Dates.format(now(),"yyyymmdd_HHMM");
  filename ="$time.jld2"
  @save filename DATA SYNC
end