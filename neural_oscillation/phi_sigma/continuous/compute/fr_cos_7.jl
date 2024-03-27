using DifferentialEquations
using Plots
using LaTeXStrings
using JLD
using Statistics
using Dates


const k_DEBUG_PRINT = false
const k_DRAW_PHASE_REALISATION = false;
const k_IS_SAVE_DATA = true;
const k_DELETE_TRANSIENT = false;
const k_DELETE_UNSTABLE = false;
const k_PRINT_ITERATION = false;

const DATA_TAKE_ERROR = 0.25;

name = "fr_cos_7"
N1 = 1
N2 = 1
const DELTA = 0.5


const NUM = 2;
const K = -500
const PAR_N = [N1, N2];
const D_MAX =  0.5
const D_ACCURACY =  0.005
const gamma1 = 1.01
const gamma2 = 1.01+DELTA
const SYNC_ERROR =  0.25;
const sigma_MAX = pi
const sigma_ACCURACY =  0.005
sigma_LIST = 0:sigma_ACCURACY:sigma_MAX
D_LIST = 0:D_ACCURACY:D_MAX
D_NUM = length(D_LIST)
sigma_NUM = length(sigma_LIST)

DATA = [zeros(4) for _ in 1:(D_NUM*sigma_NUM)]
W = [zeros(2) for _ in 1:(D_NUM*sigma_NUM)]

function eqn!(du, u, p, t)#u - это тета
  d, sigma, g, n = p
  f = g .- cos.(u ./ n)
  exch = 1 ./ (1 .+ ℯ.^(K*(cos(sigma).-sin.(u))))
  du .= f - d*exch
end

function FREQ_SYNC(DATA, gamma1, gamma2, PAR_N, D_LIST, sigma_LIST)
  Threads.@threads   for k in eachindex(sigma_LIST)
    # num_of_iterations = D_NUM*sigma_NUM
    sigma = sigma_LIST[k];
    a = 8000;
    b = 10000;

    for m in eachindex(D_LIST)    
        d = D_LIST[m]
                
        tspan = (a, b)

        p = (d, sigma, [gamma1, gamma2],  PAR_N);
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

        w_s, w_b, err = SYNC_PAIR(T, Y, PAR_N)

        ratio_s = w_s[1] / w_s[2];
        ratio_b = w_b[1] / w_b[2];

        if (err != 0)
          DATA[m] = [d, sigma,  -err, -err]
          continue;
        end

        
        DATA[m + (k-1)*D_NUM] = [d, sigma, ratio_s, ratio_b]
        W[m] = vcat(w_s, w_b)
        # it = m + (k-1)*D_NUM
        # println("Iteration $it of $num_of_iterations")
      end
    end
end

function SYNC_PAIR(T, Y, PAR_N)
  Y = reduce(vcat,transpose.(Y))
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

  global ws_1 = 2*pi./avg(Times_SP_1);
  global ws_2 = 2*pi./avg(Times_SP_2);

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

function DRAW(T, Y, gamma1, gamma2, d, PAR_N)
  Y = [mod.(y, 2 * pi) for y in Y]
  n1 = PAR_N[1];
  n2 = PAR_N[2];
  plot(T, getindex.(Y, 1), label=L"n_1 = %$n1, \gamma_{1}=%$gamma1")
  plot!(T, getindex.(Y, 2), label=L"n_2 = %$n2 , \gamma_{2}=%$gamma2")
  title!(L"d = %$d")
  ylims!(0,  2*pi)
  xlabel!(L"t")
  ylabel!(L"\varphi")
end

FREQ_SYNC(DATA, gamma1, gamma2, PAR_N, D_LIST, sigma_LIST);


if k_IS_SAVE_DATA 
  times = Dates.format(now(),"__yyyymmdd_HHMM");
  filename ="$name$times.jld2"
  @save filename DATA W
end