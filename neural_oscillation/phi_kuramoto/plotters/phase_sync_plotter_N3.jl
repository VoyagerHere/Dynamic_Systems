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
N3 = 3;

# Change it to correct number of 
# points to plot area
accuracy = 30;

@load "pi_8__3_3_3__20240411_2203.jld2" DATA SYNC DEATH

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

cmap_sp = PythonPlot.matplotlib.cm.get_cmap("summer")
cmap_brst = PythonPlot.matplotlib.cm.get_cmap("cool")
cmap_dth = PythonPlot.matplotlib.cm.get_cmap("Oranges")

cmap_sp_one = PythonPlot.matplotlib.cm.get_cmap("pink")
cmap_brst_one = PythonPlot.matplotlib.cm.get_cmap("winter")


Burst_Sync_ALL = findall(x -> x ==  1, SYNC_MATR[:,5])
Spike_Sync_ALL = findall(x -> x ==  1, SYNC_MATR[:,6])

Burst_Sync_23 = findall(x -> x ==  1, SYNC_MATR[:,3])
Spike_Sync_23 = findall(x -> x ==  1, SYNC_MATR[:,4])

Burst_Sync_12 = findall(x -> x ==  1, SYNC_MATR[:,1])
Spike_Sync_12 = findall(x -> x ==  1, SYNC_MATR[:,2])


Burst_Sync_ALL_ONE = findall(x -> x ==  2, SYNC_MATR[:,5])
Spike_Sync_ALL_ONE = findall(x -> x ==  2, SYNC_MATR[:,6])

Burst_Sync_23_ONE = findall(x -> x ==  2, SYNC_MATR[:,3])
Spike_Sync_23_ONE = findall(x -> x ==  2, SYNC_MATR[:,4])

Burst_Sync_12_ONE = findall(x -> x ==  2, SYNC_MATR[:,1])
Spike_Sync_12_ONE = findall(x -> x ==  2, SYNC_MATR[:,2])

DEAD1 = findall(x -> x ==  1, DEATH_MATR[:,1])
DEAD2 = findall(x -> x ==  1, DEATH_MATR[:,2])
DEAD3 = findall(x -> x ==  1, DEATH_MATR[:,3])
DEAD_ALL = DEAD3
QUASI_PERIODIC = findall(x -> x ==  0, DEAD_ALL)


QUASI_PERIODIC = intersect(
  findall(x -> x ==  0, SYNC_MATR[:,5]),
  findall(x -> x ==  0, SYNC_MATR[:,6]))


  if(length(QUASI_PERIODIC) > accuracy)
    scatter(D_VEC_1[QUASI_PERIODIC], DELTA_VEC[QUASI_PERIODIC], color="lavender", s=size_sc, label=L"$S$%$counter_field: Q-P")
    global counter_field+=1;
  end
  if (length(DEAD1) != length(DEAD2))
    if(length(DEAD1) > accuracy)
      scatter(D_VEC_1[DEAD1], DELTA_VEC[DEAD1], s=size_sc, color=cmap_dth(0.4), label=L"$S$%$counter_field: D $\varphi_1$")
      global counter_field+=1;
    end
  end
  
  if (length(DEAD2) != length(DEAD3))
    if(length(DEAD2) > accuracy)
      scatter(D_VEC_1[DEAD2], DELTA_VEC[DEAD2], s=size_sc, color=cmap_dth(0.6), label=L"$S$%$counter_field: D $\varphi_1, \varphi_2$")
      global counter_field+=1;
    end
  end

  
if(length(DEAD_ALL) > accuracy)
  scatter(D_VEC_1[DEAD_ALL], DELTA_VEC[DEAD_ALL], s=size_sc, color=cmap_dth(0.8), label=L"$S$%$counter_field: D $\varphi_1, \varphi_2, \varphi_3$")
  global counter_field+=1;
end

if (Burst_Sync_12 != Spike_Sync_12)
  if(length(Burst_Sync_12) > accuracy)
    burst_color = cmap_brst(0.8)
    scatter(D_VEC_1[Burst_Sync_12], DELTA_VEC[Burst_Sync_12], s=size_sc, color=burst_color, label=L"$S$%$counter_field: BS 1-2")
    global counter_field+=1;
  end
end
if(length(Spike_Sync_12) > accuracy)
  spike_color = cmap_sp(0.8)
  scatter(D_VEC_1[Spike_Sync_12], DELTA_VEC[Spike_Sync_12], s=size_sc, color=spike_color, label=L"$S$%$counter_field: SS 1-2")
  global counter_field+=1;
end

if (Burst_Sync_23 != Spike_Sync_23)
  
  if(length(Burst_Sync_23) > accuracy)
    burst_color = cmap_brst(0.4)
    scatter(D_VEC_1[Burst_Sync_23], DELTA_VEC[Burst_Sync_23], s=size_sc, color=burst_color, label=L"$S$%$counter_field: BS 2-3")
    global counter_field+=1;
  end
end
if(length(Spike_Sync_23) > accuracy)
  spike_color = cmap_sp(0.4)
  scatter(D_VEC_1[Spike_Sync_23], DELTA_VEC[Spike_Sync_23], s=size_sc, color=spike_color, label=L"$S$%$counter_field: SS 2-3")
  global counter_field+=1;
end

if (Burst_Sync_ALL != Spike_Sync_ALL)
  if(length(Burst_Sync_ALL) > accuracy)
    burst_color = cmap_brst(0)
    scatter(D_VEC_1[Burst_Sync_ALL], DELTA_VEC[Burst_Sync_ALL], s=size_sc, color=burst_color, label=L"$S$%$counter_field: GBS")
    global counter_field+=1;
  end
end
if(length(Spike_Sync_ALL) > accuracy)
  spike_color = cmap_sp(0)
  scatter(D_VEC_1[Spike_Sync_ALL], DELTA_VEC[Spike_Sync_ALL], s=size_sc, color=spike_color, label=L"$S$%$counter_field: GSS")
  global counter_field+=1;
end



if (Burst_Sync_12_ONE != Spike_Sync_12_ONE)
  if(length(Burst_Sync_12_ONE) > accuracy)
    burst_color = cmap_brst_one(0.2)
    scatter(D_VEC_1[Burst_Sync_12_ONE], DELTA_VEC[Burst_Sync_12_ONE], s=size_sc, color=burst_color, label=L"$S$%$counter_field: BS A-P 1-2")
    global counter_field+=1;
  end
end
if(length(Spike_Sync_12_ONE) > accuracy)
  spike_color = cmap_sp_one(0.2)
  scatter(D_VEC_1[Spike_Sync_12_ONE], DELTA_VEC[Spike_Sync_12_ONE], s=size_sc, color=spike_color, label=L"$S$%$counter_field: SS A-P 1-2")
  global counter_field+=1;
end

if (Burst_Sync_23_ONE != Spike_Sync_23_ONE)
  
  if(length(Burst_Sync_23_ONE) > accuracy)
    burst_color = cmap_brst_one(0.5)
    scatter(D_VEC_1[Burst_Sync_23_ONE], DELTA_VEC[Burst_Sync_23_ONE], s=size_sc, color=burst_color, label=L"$S$%$counter_field: BS A-P 2-3")
    global counter_field+=1;
  end
end
if(length(Spike_Sync_23_ONE) > accuracy)
  spike_color = cmap_sp_one(0.3)
  scatter(D_VEC_1[Spike_Sync_23_ONE], DELTA_VEC[Spike_Sync_23_ONE], s=size_sc, color=spike_color, label=L"$S$%$counter_field: SS A-P 2-3")
  global counter_field+=1;
end

if (Burst_Sync_ALL_ONE != Spike_Sync_ALL_ONE)
  if(length(Burst_Sync_ALL_ONE) > accuracy)
    burst_color = cmap_brst_one(1.0)
    scatter(D_VEC_1[Burst_Sync_ALL_ONE], DELTA_VEC[Burst_Sync_ALL_ONE], s=size_sc, color=burst_color, label=L"$S$%$counter_field: GBS A-P")
    global counter_field+=1;
  end
end
if(length(Spike_Sync_ALL_ONE) > accuracy)
  spike_color = cmap_sp_one(0.5)
  scatter(D_VEC_1[Spike_Sync_ALL_ONE], DELTA_VEC[Spike_Sync_ALL_ONE], s=size_sc, color=spike_color, label=L"$S$%$counter_field: GSS A-P")
  global counter_field+=1;
end




legend(bbox_to_anchor=(1, 1.015), loc="upper left", fontsize=12)
# legend(loc="lower right", fontsize=16, framealpha=1)
# legend(loc="upper left", fontsize=16, framealpha=1)


# for handle in lgnd.legend_handles
#   handle.set_markersize([6.0])
# end
# xlim(2.62,  2.7)
title(L"$n_1$ = %$N1, $n_2$ = %$N2, $n_3$ = %$N3, $\alpha$ = %$alpha_txt", fontsize=20)
xlabel(L"d", fontsize=20)
ylabel(L"\Delta", fontsize=20)