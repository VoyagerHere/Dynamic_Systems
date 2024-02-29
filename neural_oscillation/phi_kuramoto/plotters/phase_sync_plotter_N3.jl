# Remarks
# 1) Ratio of sync is not represented
# 2) Use only d1 in equation
#
# 3) Delta is equal for both phase 
#   phi2 = phi_1 + delta
#   phi_3 = phi_2 + delta 

using DifferentialEquations
using LaTeXStrings
using JLD
using Statistics
using Dates

using PythonPlot
pygui(true)

alpha_txt = "π/8"
N1 = 3;
N2 = 3;
N3 = 3

# Change it to correct number of 
# points to plot area
accuracy = 0;

@load "title.jld2" DATA SYNC DEATH

size = length(DATA)
DATA = reduce(vcat,transpose.(DATA))
SYNC_MATR = Int.(reduce(vcat,transpose.(SYNC)))
DEATH_MATR = Int.(reduce(vcat,transpose.(DEATH)))


D_VEC_1 = DATA[:,1]
D_VEC_2 = DATA[:,2]
RATIO_VEC_12 = DATA[:,3]
RATIO_VEC_23 = DATA[:,4]
DELTA_VEC = DATA[:,5]

counter_field = 1;

PythonPlot.matplotlib.rcParams["font.size"] = 14
size_sc = 4


Burst_Sync_ALL = findall(x -> x ==  1, SYNC_MATR[:,5])
Spike_Sync_ALL = findall(x -> x ==  1, SYNC_MATR[:,6])

Burst_Sync_23 = findall(x -> x ==  1, SYNC_MATR[:,3])
Spike_Sync_23 = findall(x -> x ==  1, SYNC_MATR[:,4])

Burst_Sync_12 = findall(x -> x ==  1, SYNC_MATR[:,1])
Spike_Sync_12 = findall(x -> x ==  1, SYNC_MATR[:,2])

DEAD1 = findall(x -> x ==  1, DEATH_MATR[:,1])
DEAD2 = findall(x -> x ==  1, DEATH_MATR[:,2])
DEAD3 = findall(x -> x ==  1, DEATH_MATR[:,3])
DEAD_ALL = intersect(DEAD1, DEAD2, DEAD3)
QUASI_PERIODIC = findall(x -> x ==  0, DEAD_ALL)


QUASI_PERIODIC = intersect(
  findall(x -> x ==  0, SYNC_MATR[:,5]),
  findall(x -> x ==  0, SYNC_MATR[:,6]))


if(length(QUASI_PERIODIC) > accuracy)
  scatter(D_VEC_1[QUASI_PERIODIC], DELTA_VEC[QUASI_PERIODIC], color = "grey", s=size_sc, label=L"$S_{counter_field}$: Q-P")
  global counter_field+=1;
end

if (Burst_Sync_12 != Spike_Sync_12)
  if(length(Burst_Sync_12) > accuracy)
    scatter(D_VEC_1[Burst_Sync_12], DELTA_VEC[Burst_Sync_12], s=size_sc, label=L"$S_{counter_field}$: BS 1-2")
    global counter_field+=1;
  end
end
if(length(Spike_Sync_12) > accuracy)
  scatter(D_VEC_1[Spike_Sync_12], DELTA_VEC[Spike_Sync_12], s=size_sc, label=L"$S_{counter_field}$: SS 1-2")
  global counter_field+=1;
end

if (Burst_Sync_23 != Spike_Sync_23)
  if(length(Burst_Sync_23) > accuracy)
    scatter(D_VEC_1[Burst_Sync_23], DELTA_VEC[Burst_Sync_23], s=size_sc, label=L"$S_{counter_field}$: BS 2-3")
    global counter_field+=1;
  end
end
if(length(Spike_Sync_23) > accuracy)
  scatter(D_VEC_1[Spike_Sync_23], DELTA_VEC[Spike_Sync_23], s=size_sc, label=L"$S_{counter_field}$: SS 2-3")
  global counter_field+=1;
end

if (Burst_Sync_ALL != Spike_Sync_ALL)
  if(length(Burst_Sync_ALL) > accuracy)
    scatter(D_VEC_1[Burst_Sync_ALL], DELTA_VEC[Burst_Sync_ALL], s=size_sc, label=L"$S_{counter_field}$: GBS")
    global counter_field+=1;
  end
end
if(length(Spike_Sync_ALL) > accuracy)
  scatter(D_VEC_1[Spike_Sync_ALL], DELTA_VEC[Spike_Sync_ALL], s=size_sc, label=L"$S_{counter_field}$: GBS")
  global counter_field+=1;
end

# legend(loc="lower right", fontsize=10)
legend(bbox_to_anchor=(1, 1.015), loc="upper left", fontsize=10)

# for handle in lgnd.legend_handles
#   handle.set_markersize([6.0])
# end

title(L"$n_1$ = %$N1, $n_2$ = %$N2, $n_3$ = %$N3, $\alpha$ = %$alpha_txt")
xlabel(L"d")
ylabel(L"\Delta")