##########################################################
## Parameter-Parallelism with GPU Ensemble Methods 
## https://docs.sciml.ai/DiffEqGPU/stable/getting_started/
## https://docs.sciml.ai/DiffEqGPU/stable/tutorials/gpu_ensemble_basic/#lorenz
###########################################################


# TODO:
# ADD delete unstable - done
# TRY delete unstable new algorithm - done (from phase realisation)
# ADD version for CUDA
# ADAPTIVE GRID - done


using DiffEqGPU, OrdinaryDiffEq, StaticArrays, CUDA
using DifferentialEquations
using Plots
using LaTeXStrings
using JLD
using Statistics

Plots.scalefontsizes()
Plots.scalefontsizes(1.5)

const k_ENABLE_ADAPTIVE_GRID = true;
const k_DEBUG_PRINT = true
const k_DRAW_PHASE_REALISATION = false;
const DATA_TAKE_ERROR = 0.05;

global a = 1000;
global b = 2000;

# For ADAPTIVE_GRID
const init_b = 2000;
const b_step = 1000;
const ADAPTIVE_SET_ERROR = 10;

const SPIKE_ERROR =  10


N1 = 1
N2 = 1
const NUM = 2;
const NUM_OF_COMPUTE_RES = 3;
global PAR_N = [N1, N2];
const D_MAX =  0.07
const D_ACCURACY =  0.0001
const G_NUM = 500
const SYNC_ERROR =  0.05
const GStart =  1.01
const DELTA_MAX_VAL =  0.025
G_LIST = range(GStart, stop=GStart + DELTA_MAX_VAL, length=G_NUM)
D_LIST = 0:D_ACCURACY:D_MAX
D_NUM = length(D_LIST)
DATA = [zeros(NUM_OF_COMPUTE_RES) for _ in 1:(D_NUM*G_NUM)]
SYNC = [zeros(NUM) for _ in 1:(D_NUM*G_NUM)]


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
    param_num = 5
    PARAM_GRID = zeros(zeros(param_num) for _ in 1:(D_NUM*G_NUM))
    G1 = GStart;
    
    # Fill param grid to send it on GPU
    
    for k in eachindex(G_LIST)
      G2 = G_LIST[k];
      for m in eachindex(D_LIST)
        d = D_LIST[m]
        PARAM_GRID[m + (k-1)*D_NUM] = [d, ALPHA, [G1, G2], PAR_N, NUM];
      end
    end

    y0 = @SVector[0; 0]
    tspan = (a, b)
    p = @SVector[d, ALPHA, [G1, G2], PAR_N, NUM];
    prob = ODEProblem(eqn!, y0, tspan, p)
    prob_func = (prob, i, repeat) -> remake(prob, p = (@SVector PARAM_GRID[i]))
    monteprob = EnsembleProblem(prob, prob_func = prob_func, safetycopy = false)
    sol = solve(monteprob, GPUTsit5(), EnsembleGPUKernel(CUDA.CUDABackend()),
    trajectories = D_NUM*G_NUM,
    saveat = 1.0f0)

    for k in eachindex(G_LIST)
      for m in eachindex(D_LIST)
        Y = sol.u[m + (k-1)*D_NUM];
        T = sol.t[m + (k-1)*D_NUM];
        DIFF_SP, DIFF_BS, ratio, err = SYNC_PAIR(T, Y, PAR_N, SPIKE_ERROR)
        sync = [0, 0]
        if (err == 0)
          sync[1] = IS_SYNC(DIFF_BS, SYNC_ERROR);
          sync[2] = IS_SYNC(DIFF_SP, SYNC_ERROR);
        end
        delta = G2 - G1
        # NUM_OF_COMPUTE_RES
        DATA[m + (k-1)*D_NUM] = [d, ratio, delta]
        SYNC[m + (k-1)*D_NUM] = sync;
      end
    end
end

function SYNC_PAIR(T, Y, PAR_N, error)
  Y = reduce(vcat,transpose.(Y))
  SPIKES1, err1 = FIND_SPIKES(Y[:,1], PAR_N[1])
  SPIKES2, err2 = FIND_SPIKES(Y[:,2], PAR_N[2])

  unstbl_1 = DELETE_UNSTBL(SPIKES1, err1, PAR_N[1], error)
  unstbl_2 = DELETE_UNSTBL(SPIKES2, err2, PAR_N[2], error)

  SPIKES1 = SPIKES1[unstbl_1:end];
  SPIKES2 = SPIKES2[unstbl_2:end];

  k_DEBUG_PRINT && println("Spikes in 1: ", length(SPIKES1))
  k_DEBUG_PRINT && println("Spikes in 2: ", length(SPIKES2))

  err = handle_errors(Bool(err1), Bool(err2));
  if err != 0
      DIFF_SP = 0
      DIFF_BS = 0
      ratio = -err
      return (DIFF_SP, DIFF_BS, ratio, err)
  end

  # UNSTBL1 = DELETE_UNSTBL(A1, err1, n, 1, error)
  # A1 = A1[UNSTBL1:end]
  BURSTS1 = FIND_BURST(SPIKES1, PAR_N[1])
  BURSTS2 = FIND_BURST(SPIKES2, PAR_N[2])

  k_ENABLE_ADAPTIVE_GRID && (ADAPTIVE_GRID(minimum([length(BURSTS1), length(BURSTS2)]), false))

  
  ratio = FIND_RATIO(BURSTS1, BURSTS2)

  DIFF_SP = FIND_DIFF(SPIKES1, SPIKES2, T)
  DIFF_BS = FIND_DIFF(BURSTS1, BURSTS2, T)
  return (DIFF_SP, DIFF_BS, ratio, err)
end

# function ADAPTIVE_GRID(num_of_bursts, reset)
#   global b;
#   if (reset == true)
#     b = init_b;
#   else   
#     if (num_of_bursts < ADAPTIVE_SET_ERROR)
#       b += b_step;
#     end
#   end
# end


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
  Y = Y .% 2*pi;
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

function FIND_NEAR_POINTS(POINTS)
  i = 1;
    while i < length(POINTS)
        if POINTS[i+1] - POINTS[i] ==  1
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

function DETECT_DEATH(errors)
    if errors[1] && errors[2]
        err = 11
    elseif errors[2]
        err = 10
    elseif errors[1]
        err = 01
    else
        err = 0
    end
    return err
end

function DRAW(T, Y, G, D, PAR_N)
    Y = [mod.(y, 2 * pi) for y in Y]
    for i in eachindex(PAR_N)
      n_cur = PAR_N[i];
      g_cur = G[i];
      plot(T, getindex.(Y, i), label=L"n_%$i = %$n_cur , \gamma_{%$i}=%$g_cur")
    end
    title!(L"%$D")
    ylims!(0,  2*pi)
    xlabel!(L"t")
    ylabel!(L"\varphi")
end




DATA = PHASE_SYNC(DATA, SYNC, GStart, PAR_N, NUM, G_LIST, D_LIST, SPIKE_ERROR, ALPHA);


# save_object("DATA.jld", DATA)
# save_object("SYNC.jld", SYNC)