using DifferentialEquations
using Plots
using LaTeXStrings
using JLD
using Statistics
using Dates

plotlyjs()

Plots.scalefontsizes()
Plots.scalefontsizes(1.5)

G1 = 1.01;
alpha_txt = "π/8"
N1 = 2;
N2 = 2;

@load "pi_8_2.jld2" DATA SYNC

size = length(DATA)
DATA = reduce(vcat,transpose.(DATA))
SYNC = Int.(reduce(vcat,transpose.(SYNC)))

D_VEC = DATA[:,1]
RATIO_VEC = Int.(DATA[:,2])
DELTA_VEC = DATA[:,3]

unique_ratios = unique(RATIO_VEC)
Burst_Sync = findall(x -> x ==  1, SYNC[:,1])
Spike_Sync = findall(x -> x ==  1, SYNC[:,2])

counter_field = 1;
for ratio in unique_ratios
  if (ratio == NaN)
    continue
  end
  ratio_ind = findall(x -> x == ratio, RATIO_VEC)

  if (ratio == 0)
    scatter!(D_VEC[ratio_ind], DELTA_VEC[ratio_ind], label=L"S_%$counter_field: Quasi-Periodic")
    global counter_field+=1;
  elseif (ratio == -1)
    scatter!(D_VEC[ratio_ind], DELTA_VEC[ratio_ind], label=L"S_%$counter_field: Death \phi_1")
    global counter_field+=1;
  elseif (ratio == -2)
    scatter!(D_VEC[ratio_ind], DELTA_VEC[ratio_ind], label=L"S_%$counter_field: Death \phi_2")
    global counter_field +=1;
  elseif (ratio == -3)
    scatter!(D_VEC[ratio_ind], DELTA_VEC[ratio_ind], label=L"S_%$counter_field: Death \phi_1, \phi_2")
    global counter_field +=1;
  else
    Burst_p = intersect(Burst_Sync, ratio_ind)
    Spike_p = intersect(Spike_Sync, ratio_ind)
    scatter!(D_VEC[Burst_p], DELTA_VEC[Burst_p], label=L"S_%$counter_field Burst sync  1:%$ratio") 
    global counter_field+=1
    scatter!(D_VEC[Spike_p], DELTA_VEC[Spike_p], label=L"S_%$counter_field Spike sync  1:%$ratio")
    global counter_field+=1
  end
end
title!(L"n_1 = %$N1, n_2 = %$N2,  \alpha = %$alpha_txt, \gamma_{1}=1.01, \gamma_{2}=\gamma_{1}+\DELTA_VEC")
xlabel!(L"d")
ylabel!(L"\delta")