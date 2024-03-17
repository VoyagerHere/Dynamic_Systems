using DifferentialEquations
using Plots
using LaTeXStrings
using JLD
using Statistics
using Dates


const k_ENABLE_ADAPTIVE_GRID = true;
const k_DEBUG_PRINT = false
const k_DRAW_PHASE_REALISATION = false;
const k_IS_SAVE_DATA = true;
const k_DELETE_TRANSIENT = false; 
const k_DELETE_UNSTABLE = true;
const k_PRINT_ITERATION = false;

const DATA_TAKE_ERROR = 0.25;
const b_init = 14000

# For ADAPTIVE_GRID
const b_step = 3000;
const ADAPTIVE_SET_ERROR = 10;

# For DELETE_UNSTABLE
const SPIKE_ERROR =  20

name = "pi_8__2_2_2"
N1 = 2;
N2 = 2;
N3 = 2;
const ALPHA = pi/8

const NUM = 3;
const PAR_N = [N1, N2, N3];
const D_MAX =  0.05
const D_ACCURACY =  0.0001
const G_NUM = 600
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
    G1 = GStart;
    death_state = [0, 0, 0]

    Threads.@threads for k in eachindex(G_LIST)
      G2 = G_LIST[k]
      G3 = G2 + DELTA;
      a = 10000;
      b = b_init;
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

        DIFF_SP_12, DIFF_BS_12, ratio_12, err_12, b1, len_1_2 = SYNC_PAIR(T, Y, PAR_N, SPIKE_ERROR, 1, 2, b)
        DIFF_SP_23, DIFF_BS_23, ratio_23, err_23, b2, len_2_3 = SYNC_PAIR(T, Y, PAR_N, SPIKE_ERROR, 2, 3, b)
        b = max(b1, b2);
        err = vcat(err_12, err_23);
        deleteat!(err, 2);

        len = vcat(len_1_2, len_2_3);
        deleteat!(len, 2);

        sync = zeros(NUM * 2);
        delta = G2 - G1

        if ((sum(err) > 0))
          for i in eachindex(err)
            if((err[i] == 1) && (death_state[i] == 0))
              b = b_init;
              death_state[i] = 1;
              k_DEBUG_PRINT && println("Neuron - ", i, " DEAD")
            end
          end
        end

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

        @inbounds DATA[m + (k-1)*D_NUM] = [d1, d2, ratio_12, ratio_23, delta]
        @inbounds SYNC[m + (k-1)*D_NUM] = sync;
        @inbounds DEATH[m + (k-1)*D_NUM] = err;
      end
    end
end

function SYNC_PAIR(T, Y, PAR_N, error, ind1, ind2, b)
  Y = reduce(vcat,transpose.(Y))
  SPIKES1, err1 = FIND_SPIKES(Y[:,ind1], PAR_N[ind1])
  SPIKES2, err2 = FIND_SPIKES(Y[:,ind2], PAR_N[ind2])

  if (k_DELETE_UNSTABLE)
    unstbl_1 = DELETE_UNSTBL(SPIKES1, err1, PAR_N[1], error)
    unstbl_2 = DELETE_UNSTBL(SPIKES2, err2, PAR_N[2], error)

    SPIKES1 = SPIKES1[unstbl_1:end];
    SPIKES2 = SPIKES2[unstbl_2:end];
  end

  err = [err1, err2];

  if sum(err) > 0 
      DIFF_SP = 0
      DIFF_BS = 0
      ratio = -1
      return (DIFF_SP, DIFF_BS, ratio, err, b, [div(length(SPIKES1), PAR_N[ind1]),div(length(SPIKES2), PAR_N[ind2])])
  end

  BURSTS1 = FIND_BURST(SPIKES1, PAR_N[ind1])
  BURSTS2 = FIND_BURST(SPIKES2, PAR_N[ind2])

  k_DEBUG_PRINT && println("Bursts in 1: ", length(BURSTS1))
  k_DEBUG_PRINT && println("Bursts in 2: ", length(BURSTS2))

  len = [length(BURSTS1), length(BURSTS2)]

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
  if err == 1
      return UNSTBL
  end
  B = zeros(Int64, div(length(SPIKES), n))
  for m in n:n:length(SPIKES)
    B[div(m, n)] = SPIKES[m]
  end
  if length(B) < 2*unstable_err
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
      rel_change_3 = abs((Y[i][3] - Y[i-1][3]) / Y[i-1][3])
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

function DRAW(T, Y, G1, G2, G3, D, PAR_N)
  Y = [mod.(y, 2 * pi) for y in Y]
  n1 = PAR_N[1];
  n2 = PAR_N[2];
  n3 = PAR_N[3];
  plot(T, getindex.(Y, 1), label=L"n_1 = %$n1, \gamma_{1}=%$G1")
  plot!(T, getindex.(Y, 2), label=L"n_2 = %$n2, \gamma_{2}=%$G2")
  plot!(T, getindex.(Y, 3), label=L"n_3 = %$n3, \gamma_{2}=%$G3")
  title!(L"d = %$D")
  ylims!(0,  2*pi)
  xlabel!(L"t")
  ylabel!(L"\varphi")
end

@time PHASE_SYNC(DATA, SYNC, GStart, PAR_N, NUM, G_LIST, D_LIST, SPIKE_ERROR, ALPHA);

if k_IS_SAVE_DATA 
  times = Dates.format(now(),"__yyyymmdd_HHMM");
  filename ="$name$times.jld2"
  @save filename DATA SYNC DEATH
end