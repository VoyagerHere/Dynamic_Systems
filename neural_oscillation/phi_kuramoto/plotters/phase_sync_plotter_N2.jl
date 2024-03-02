using DifferentialEquations
using LaTeXStrings
using JLD
using Statistics
using Dates

using PythonPlot
pygui(true)

name = "untitiled"
alpha_txt = "π/8"
N1 = 3;
N2 = 3;

@load "$name.jld2" DATA SYNC

size = length(DATA)
DATA = reduce(vcat,transpose.(DATA))
SYNC = Int.(reduce(vcat,transpose.(SYNC)))

D_VEC = DATA[:,1]
RATIO_VEC = Int.(DATA[:,2])
DELTA_VEC = DATA[:,3]

Burst_Sync = findall(x -> x ==  1, SYNC[:,1])
Spike_Sync = findall(x -> x ==  1, SYNC[:,2])

unique_ratios = unique(RATIO_VEC)
unique_ratios = sort(unique_ratios, rev=true)
counter_field = 1;

accuracy = 10;
PythonPlot.matplotlib.rcParams["font.size"] = 14
size_sc = 4

unique_ratios = sort!(unique_ratios, alg=InsertionSort);

for ratio in unique_ratios
  if (ratio == NaN)
    continue
  end 
  ratio_ind = findall(x -> x == ratio, RATIO_VEC)
  if (ratio == 0)
    scatter(D_VEC[ratio_ind], DELTA_VEC[ratio_ind], s=size_sc, color = "gray", label=L"$S_{%$counter_field}$: Q-P")
    global counter_field+=1;
  elseif (ratio == -1)
    scatter(D_VEC[ratio_ind], DELTA_VEC[ratio_ind], s=size_sc, color = "darkgrey", label=L"$S_{%$counter_field}$: D $\varphi_1$")
    global counter_field += 1;
  elseif (ratio == -2)
    scatter(D_VEC[ratio_ind], DELTA_VEC[ratio_ind], s=size_sc, color = "slategrey", label=L"$S_{%$counter_field}$: D $\varphi_2$")
    global counter_field +=1;
  elseif (ratio == -3)
    scatter(D_VEC[ratio_ind], DELTA_VEC[ratio_ind], s=size_sc, color = "black", label=L"$S_{%$counter_field}$: D $\varphi_1, \varphi_2$")
    global counter_field +=1;
  else
    Burst_p = intersect(Burst_Sync, ratio_ind)
    Spike_p = intersect(Spike_Sync, ratio_ind)
    if (length(Spike_p) > accuracy) && (length(Burst_p) > accuracy)
      if (Spike_p == Burst_p) 
        scatter(D_VEC[Spike_p], DELTA_VEC[Spike_p], s=size_sc, label=L"$S_{%$counter_field}$: SS 1:%$ratio")
        global counter_field+=1
      else
        scatter(D_VEC[Burst_p], DELTA_VEC[Burst_p], s=size_sc, label=L"$S_{%$counter_field}$: BS  1:%$ratio") 
        global counter_field+=1
        scatter(D_VEC[Spike_p], DELTA_VEC[Spike_p], s=size_sc, label=L"$S_{%$counter_field}$: SS 1:%$ratio")
        global counter_field+=1
      end
    end
  end
end
# legend(loc="lower right", fontsize=10)
legend(bbox_to_anchor=(1, 1.015), loc="upper left", fontsize=10)

title(L"$n_1$ = %$N1, $n_2$ = %$N2,  $\alpha$ = %$alpha_txt")
xlabel(L"d")
ylabel(L"\Delta")