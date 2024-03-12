using DifferentialEquations
using Plots
using LaTeXStrings
using JLD
using Statistics
using Dates



const k_ENABLE_ADAPTIVE_GRID = true;
const k_DRAW_PHASE_REALISATION = false;
const k_IS_SAVE_DATA = true;
const k_DELETE_TRANSIENT = false;
const k_DELETE_UNSTABLE = false;
const k_ADAPTIVE_SOL_POINTS = true;


const k_DEBUG_PRINT = true
const k_PRINT_ITERATION = false;
const k_PRINT_ERROR = true;

const DATA_TAKE_ERROR = 0.25;

const ADAPTIVE_SET_ERROR = 10;
const SPIKE_ERROR =  0

name = "pi_2_3__3_3_3"
const ALPHA = 2*pi / 3
N1 = 3
N2 = 3
N3 = 3


const b_init = 9000;
const b_step = 1000;


const NUM = 3;
global PAR_N = [N1, N2, N3];
const D_MAX =  0.05
const D_ACCURACY =  0.00001
const SYNC_ERROR =  0.05
const GStart =  1.01
const DELTA =  0.025


D_LIST = 0:D_ACCURACY:D_MAX
D_NUM = length(D_LIST)

G1 = 1.01;
G2 = G1 + DELTA;
G3 = G2 + DELTA;


const NUM_OF_COMPUTE_RES = 6;
DATA = [zeros(NUM_OF_COMPUTE_RES) for _ in 1:D_NUM]
W = [zeros(4) for _ in 1:D_NUM]
DEATH = [zeros(NUM) for _ in 1:(D_NUM)]



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
  death_state = [0, 0, 0]

  a = 8000;
  b = b_init;

  for m in eachindex(D_LIST)
    global d1 = D_LIST[m]
    global d2 = d1;    
    tspan = (a, b)
    
    p = ([d1, d2], ALPHA, [G1, G2, G3], PAR_N, NUM);
    y0 = [0; 0; 0]

    prob = ODEProblem(eqn!, y0, tspan, p)
    sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
    global Y = sol.u;
    global T = sol.t;

    if (k_DELETE_TRANSIENT)
      index = DELETE_TRANSIENT(Y)
      start = T[index];
      y0 = [start, start, start]

      prob = ODEProblem(eqn!, y0, tspan, p)
      sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
      Y = sol.u;
      T = sol.t;
    end

    w_s_1_2, w_b_1_2, err_1, len_1_2, b_1 = SYNC_PAIR(T, Y, PAR_N, SPIKE_ERROR, 1, b)
    w_s_2_3, w_b_2_3, err_2, len_2_3, b_2 = SYNC_PAIR(T, Y, PAR_N, SPIKE_ERROR, 2, b)
    b = max(b_1, b_2);

    len = vcat(len_1_2, len_2_3);
    deleteat!(len, 2);
    global err = vcat(err_1, err_2)
    deleteat!(err, 2)

    if ((sum(err) > 0))
      for i in eachindex(err)
        if((err[i] == 1) && (death_state[i] == 0))
          b = b_init;
          death_state[i] = 1;
          k_DEBUG_PRINT && println("Neuron - ", i, " DEAD")
        end
      end
    end


    
    if (k_ADAPTIVE_SOL_POINTS && (sum(err) != 3))
      k_ADAPTIVE_TOL = false
      accuracy = 0.01;
      step = 500;
      coef = 3/2;
      while (((len[1] % PAR_N[1]) != 0) || ((len[2] % PAR_N[2] != 0))|| ((len[3] % PAR_N[3] != 0)))
        if (k_ADAPTIVE_TOL)
          saveat_points = a:accuracy:b+step
          prob = ODEProblem(eqn!, y0, (a, b+step), p)
          sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12, saveat=saveat_points)
        else
          prob = ODEProblem(eqn!, y0, (a, b+step), p)
          sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
        end
        global Y = sol.u;
        global T = sol.t;

        w_s_1_2, w_b_1_2, err_1, len_1_2, b_1 = SYNC_PAIR(T, Y, PAR_N, SPIKE_ERROR, 1, b)
        w_s_2_3, w_b_2_3, err_2, len_2_3, b_2 = SYNC_PAIR(T, Y, PAR_N, SPIKE_ERROR, 2, b)
        b = max(b_1, b_2);
    
        len = vcat(len_1_2, len_2_3);
        deleteat!(len, 2);
        err = vcat(err_1, err_2)
        deleteat!(err, 2)

        accuracy = accuracy / coef;
        step = step / coef;          
        if (accuracy < 1e-3)
          k_PRINT_ERROR && println("ERROR: Too big tolerance in D1:$d1, D2:$d2, G1:$G1, G2:$G2, len:$len")
          break;
        end
      end
    end
    

    if (sum(err) == 3)
      @inbounds DATA[m] = [d1, d2, -1, -1, -1, -1]
      @inbounds W[m] = vcat(0, 0, 0, 0)
      @inbounds DEATH[m] = err;
      continue;
    end

    ratio_s_1_2 = w_s_1_2[1] / w_s_1_2[2];
    ratio_s_2_3 = w_s_2_3[1] / w_s_2_3[2];
    ratio_b_1_2 = w_b_1_2[1] / w_b_1_2[2];
    ratio_b_2_3 = w_b_2_3[1] / w_b_2_3[2];

    @inbounds DATA[m] = [d1, d2, ratio_s_1_2, ratio_b_1_2, ratio_s_2_3, ratio_b_2_3]
    @inbounds W[m] = vcat(w_s_1_2, w_b_1_2, w_s_2_3, w_b_2_3)
    @inbounds DEATH[m] = err;
    k_PRINT_ITERATION && println("Iteration $m of $num_of_iterations")
  end
end

function SYNC_PAIR(T, Y, PAR_N, error, ind, b)
  Y = reduce(vcat,transpose.(Y))
  global SPIKES1, err1 = FIND_SPIKES(Y[:,ind], PAR_N[ind])
  global SPIKES2, err2 = FIND_SPIKES(Y[:,ind+1], PAR_N[ind+1])

  len = [length(SPIKES1), length(SPIKES2)]

  if (k_DELETE_UNSTABLE)
    unstbl_1 = DELETE_UNSTBL(SPIKES1, err1, PAR_N[ind], error)
    unstbl_2 = DELETE_UNSTBL(SPIKES1, err2, PAR_N[ind+1], error)

    untsbl_ind = max(unstbl_1, unstbl_2)

    SPIKES1 = SPIKES1[untsbl_ind:end];
    SPIKES2 = SPIKES2[untsbl_ind:end];
  end

  k_DEBUG_PRINT && println("Spikes in 1: ", length(SPIKES1))
  k_DEBUG_PRINT && println("Spikes in 2: ", length(SPIKES2))

  err = [err1, err2];

  if sum(err) > 0 
    return ([0, 0], [0, 0], err, len, b)
  end

  k_ENABLE_ADAPTIVE_GRID && (b = ADAPTIVE_GRID(minimum([length(SPIKES1), length(SPIKES1)]), b))

  Times_1 = FIND_TIMES(SPIKES1, T, PAR_N);
  Times_2 = FIND_TIMES(SPIKES2, T, PAR_N);

  Times_SP_1 = DEL_MIDDLE_BURST_INTRVL(PAR_N[ind], Times_1)
  Times_SP_2 = DEL_MIDDLE_BURST_INTRVL(PAR_N[ind], Times_2)

  avg(x) = (ones(length(x)) / length(x))'*x

  wb_1 = 2*pi./sum(Times_1);
  wb_2 = 2*pi./sum(Times_2);

  ws_1 = 2*pi./avg(Times_SP_1);
  ws_2 = 2*pi./avg(Times_SP_2);

  return ([ws_1, ws_2], [wb_1, wb_2], err, len, b)
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

function ADAPTIVE_GRID(num_of_bursts, b)
  if (num_of_bursts < ADAPTIVE_SET_ERROR)
    b = b +  b_step;
  end
  return b;
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
FREQ_SYNC(DATA, G1, G2, PAR_N, NUM, D_LIST, SPIKE_ERROR, ALPHA);


if k_IS_SAVE_DATA 
  times = Dates.format(now(),"__yyyymmdd_HHMM");
  filename ="$name$times.jld2"
  @save filename DATA W
end