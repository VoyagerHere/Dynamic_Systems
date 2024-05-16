using DifferentialEquations
using LaTeXStrings
using JLD
using Statistics
using Dates

using PythonPlot
pygui(true)

N1 = 1;
N2 = 1;

@load "ph_discr__20240417_2310.jld2" DATA SYNC

size = length(DATA)
DATA = reduce(vcat,transpose.(DATA))
SYNC = Int.(reduce(vcat,transpose.(SYNC)))

D_VEC = DATA[:,1]
RATIO_VEC = Int.(DATA[:,2])
SIGMA_VEC = DATA[:,3]

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
    scatter(SIGMA_VEC[ratio_ind], D_VEC[ratio_ind], s=size_sc, color = "lightcoral", label=L"$S$%$counter_field: Q-P")
    global counter_field+=1;
  elseif (ratio == -1)
    scatter(SIGMA_VEC[ratio_ind], D_VEC[ratio_ind], s=size_sc, color = "indianred", label=L"$S$%$counter_field: D $\varphi_1$")
    global counter_field += 1;
  elseif (ratio == -2)
    scatter(SIGMA_VEC[ratio_ind], D_VEC[ratio_ind], s=size_sc, color = "firebrick", label=L"$S$%$counter_field: D $\varphi_2$")
    global counter_field +=1;
  elseif (ratio == -3)
    scatter(SIGMA_VEC[ratio_ind], D_VEC[ratio_ind], s=size_sc, color = "darkred", label=L"$S$%$counter_field: D $\varphi_1, \varphi_2$")
    global counter_field +=1;
  else
    Burst_p = intersect(Burst_Sync, ratio_ind)
    Spike_p = intersect(Spike_Sync, ratio_ind)
    if (length(Spike_p) > accuracy) && (length(Burst_p) > accuracy)
      if (Spike_p == Burst_p) 
        scatter(SIGMA_VEC[Spike_p], D_VEC[Spike_p], s=size_sc, label=L"$S$%$counter_field: SS 1:%$ratio")
        global counter_field+=1
      else
        scatter(SIGMA_VEC[Burst_p], D_VEC[Burst_p], s=size_sc, label=L"$S$%$counter_field: BS 1:%$ratio") 
        global counter_field+=1
        scatter(SIGMA_VEC[Spike_p], D_VEC[Spike_p], s=size_sc, label=L"$S$%$counter_field: SS 1:%$ratio")
        global counter_field+=1
      end
    end
  end
end
# legend(loc="lower right", fontsize=10)
legend(bbox_to_anchor=(1, 1.015), loc="upper left", fontsize=10)

title(L"$n_1$ = %$N1, $n_2$ = %$N2")
xlabel(L"\Sigma")
ylabel(L"\d")