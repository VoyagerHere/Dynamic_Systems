using DifferentialEquations
using LaTeXStrings
using JLD
using Statistics
using Dates

using PythonPlot
pygui(true)

G1 = 0;
N1 = 0;
N2 = 0;

@load "name.jld2" DATA SYNC

size = length(DATA)
DATA = reduce(vcat,transpose.(DATA))
SYNC = Int.(reduce(vcat,transpose.(SYNC)))

D_VEC = DATA[:,1]
RATIO_VEC = Int.(DATA[:,2])
SIGMA_VEC = DATA[:,3]

unique_ratios = unique(RATIO_VEC)
Burst_Sync = findall(x -> x ==  1, SYNC[:,1])
Spike_Sync = findall(x -> x ==  1, SYNC[:,2])

counter_field = 1;

PythonPlot.matplotlib.rcParams["font.size"] = 14
size_sc = 4

unique_ratios = sort!(unique_ratios, alg=InsertionSort);

for ratio in unique_ratios
  if (ratio == NaN)
    continue
  end 
  ratio_ind = findall(x -> x == ratio, RATIO_VEC)
  if (ratio == 0)
    scatter(SIGMA_VEC[ratio_ind], D_VEC[ratio_ind], s=size_sc, color = "gray", label=L"$S_{1}$: Q-P")
    global counter_field+=1;
  elseif (ratio == -1)
    scatter(SIGMA_VEC[ratio_ind], D_VEC[ratio_ind], s=size_sc, color = "darkgrey", label=L"$S_{2}$: D $\varphi_1$")
    global counter_field += 1;
  elseif (ratio == -2)
    scatter(SIGMA_VEC[ratio_ind], D_VEC[ratio_ind], s=size_sc, color = "slategrey", label=L"$S_{3}$: D $\varphi_2$")
    global counter_field +=1;
  elseif (ratio == -3)
    scatter(SIGMA_VEC[ratio_ind], D_VEC[ratio_ind], s=size_sc, color = "black", label=L"$S_{4}$: D $\varphi_1, \varphi_2$")
    global counter_field +=1;
  else
    Burst_p = intersect(Burst_Sync, ratio_ind)
    Spike_p = intersect(Spike_Sync, ratio_ind)
    scatter(SIGMA_VEC[Burst_p], D_VEC[Burst_p], s=size_sc, label=L"$S_{%$counter_field}$: BS  1:%$ratio") 
    global counter_field+=1
    scatter(SIGMA_VEC[Spike_p], D_VEC[Spike_p], s=size_sc, label=L"$S_{%$counter_field}$: SS 1:%$ratio")
    global counter_field+=1
  end
end
# legend(loc="lower right", fontsize=10)
legend(bbox_to_anchor=(1, 1.015), loc="upper left", fontsize=10)

for handle in lgnd.legend_handles
  handle.set_markersize([6.0])
end

title(L"$n_1$ = %$N1, $n_2$ = %$N2")
xlabel(L"\Sigma")
ylabel(L"\d")