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

G1 = 0;
alpha_txt = "π/8"
N1 = 0;
N2 = 0;

# Change it to correct number of 
# points to plot area
accuracy = 0;

@load "name.jld2" DATA SYNC

size = length(DATA)
DATA = reduce(vcat,transpose.(DATA))
SYNC_MATR = Int.(reduce(vcat,transpose.(SYNC)))

D_VEC_1 = DATA[:,1]
D_VEC_2 = DATA[:,2]
RATIO_VEC_12 = DATA[:,3]
RATIO_VEC_23 = DATA[:,4]
DELTA_VEC = DATA[:,5]

unique_ratios = unique(RATIO_VEC)


counter_field = 1;

PythonPlot.matplotlib.rcParams["font.size"] = 14
size_sc = 4


Burst_Sync_ALL = findall(x -> x ==  1, SYNC[:,5])
Spike_Sync_ALL = findall(x -> x ==  1, SYNC[:,6])

Burst_Sync_23 = findall(x -> x ==  1, SYNC[:,3])
Spike_Sync_23 = findall(x -> x ==  1, SYNC[:,4])

Burst_Sync_12 = findall(x -> x ==  1, SYNC[:,1])
Spike_Sync_12 = findall(x -> x ==  1, SYNC[:,2])


global counter_field+=1;

if(length(Burst_Sync_12) > accuracy)
  scatter(D_VEC_1[Burst_Sync_12], DELTA_VEC[Burst_Sync_12], s=size_sc, label=L"$S_{1}$: BS 1-2")
  global counter_field+=1;
end

if(length(Spike_Sync_12) > accuracy)
  scatter(D_VEC_1[Spike_Sync_12], DELTA_VEC[Spike_Sync_12], s=size_sc, label=L"$S_{1}$: SS 1-2")
  global counter_field+=1;
end

if(length(Burst_Sync_23) > accuracy)
  scatter(D_VEC_1[Burst_Sync_23], DELTA_VEC[Burst_Sync_23], s=size_sc, label=L"$S_{1}$: BS 2-3")
  global counter_field+=1;
end

if(length(Spike_Sync_23) > accuracy)
  scatter(D_VEC_1[Spike_Sync_23], DELTA_VEC[Spike_Sync_23], s=size_sc, label=L"$S_{1}$: SS 2-3")
  global counter_field+=1;
end

if(length(Burst_Sync_ALL) > accuracy)
  scatter(D_VEC_1[Burst_Sync_ALL], DELTA_VEC[Burst_Sync_ALL], s=size_sc, label=L"$S_{1}$: GBS")
  global counter_field+=1;
end

if(length(Spike_Sync_ALL) > accuracy)
  scatter(D_VEC_1[Spike_Sync_ALL], DELTA_VEC[Spike_Sync_ALL], s=size_sc, label=L"$S_{1}$: GBS")
  global counter_field+=1;
end

# legend(loc="lower right", fontsize=10)
legend(bbox_to_anchor=(1, 1.015), loc="upper left", fontsize=10)

for handle in lgnd.legend_handles
  handle.set_markersize([6.0])
end

title(L"$n_1$ = %$N1, $n_2$ = %$N2,  $\alpha$ = %$alpha_txt")
xlabel(L"d")
ylabel(L"\Delta")