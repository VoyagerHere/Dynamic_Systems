using DifferentialEquations
using LaTeXStrings
using JLD
using Statistics
using Dates

using PythonPlot
pygui(true)

name = "pi_8__3_3__20240318_1441"
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

cmap_org = PythonPlot.matplotlib.cm.get_cmap("Oranges")
cmap_bl = PythonPlot.matplotlib.cm.get_cmap("Blues")
cmap_dth = PythonPlot.matplotlib.cm.get_cmap("Reds")


unique_ratios = sort!(unique_ratios, alg=InsertionSort);

for ratio in unique_ratios
  if (ratio == NaN)
    continue
  end 
  ratio_ind = findall(x -> x == ratio, RATIO_VEC)
  if (ratio == 0)
    scatter(D_VEC[ratio_ind], DELTA_VEC[ratio_ind], s=size_sc, c=cmap_dth(0.4), label=L"$S$%$counter_field: Q-P")
    global counter_field+=1;
  elseif (ratio == -1)
    scatter(D_VEC[ratio_ind], DELTA_VEC[ratio_ind], s=size_sc, c=cmap_dth(0.6), label=L"$S$%$counter_field: D $\varphi_1$")
    global counter_field += 1;
  elseif (ratio == -2)
    scatter(D_VEC[ratio_ind], DELTA_VEC[ratio_ind], s=size_sc, c=cmap_dth(0.8), label=L"$S$%$counter_field: D $\varphi_2$")
    global counter_field +=1;
  elseif (ratio == -3)
    scatter(D_VEC[ratio_ind], DELTA_VEC[ratio_ind], s=size_sc, c=cmap_dth(1), label=L"$S$%$counter_field: D $\varphi_1, \varphi_2$")
    global counter_field +=1;
  else
    Burst_p = intersect(Burst_Sync, ratio_ind)
    Spike_p = intersect(Spike_Sync, ratio_ind)
    if (length(Spike_p) > accuracy) || (length(Burst_p) > accuracy)
      if (Spike_p == Burst_p) 
        color = cmap_org(length(Spike_p) / size) # Normalize the length to a value between 0 and 1
        
        scatter(D_VEC[Spike_p], DELTA_VEC[Spike_p], s=size_sc, c=color,  label=L"$S$%$counter_field: SS 1:%$ratio")
        global counter_field+=1
      else
        if ((length(Burst_p) > accuracy))
          burst_color = cmap_bl(length(Burst_p) / size) # Normalize the length to a value between 0 and 1

          scatter(D_VEC[Burst_p], DELTA_VEC[Burst_p], s=size_sc, c=burst_color, label=L"$S$%$counter_field: BS 1:%$ratio") 
          global counter_field+=1
        end
        if ((length(Spike_p) > accuracy))
          spike_color = cmap_org(length(Spike_p) / size) # Normalize the length to a value between 0 and 1

          scatter(D_VEC[Spike_p], DELTA_VEC[Spike_p], s=size_sc, c=spike_color,label=L"$S$%$counter_field: SS 1:%$ratio")
          global counter_field+=1
        end
      end
    end
  end
end
legend(loc="lower right", fontsize=16, framealpha=1)
# legend(bbox_to_anchor=(1, 1.015), loc="upper left", fontsize=12)

title(L"$n_1$ = %$N1, $n_2$ = %$N2,  $\alpha$ = %$alpha_txt", fontsize=20)
xlabel(L"d", fontsize=20)
ylabel(L"\Delta", fontsize=20)