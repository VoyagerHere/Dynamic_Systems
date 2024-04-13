using DifferentialEquations
using LaTeXStrings
using JLD
using Statistics
using Dates

using PythonPlot
pygui(true)

name = "pi_8__3_3__20240411_2257"
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

Burst_Sync_ONE = findall(x -> x ==  2, SYNC[:,1])
Spike_Sync_ONE = findall(x -> x ==  2, SYNC[:,2])


unique_ratios = unique(RATIO_VEC)
unique_ratios = sort(unique_ratios, rev=true)
counter_field = 1;

accuracy = 40;
PythonPlot.matplotlib.rcParams["font.size"] = 14
size_sc = 4

cmap_sp = PythonPlot.matplotlib.cm.get_cmap("summer")
cmap_brst = PythonPlot.matplotlib.cm.get_cmap("cool")
cmap_dth = PythonPlot.matplotlib.cm.get_cmap("Oranges")


unique_ratios = sort!(unique_ratios, alg=InsertionSort);

for ratio in unique_ratios
  if (ratio == NaN)
    continue
  end 
  ratio_ind = findall(x -> x == ratio, RATIO_VEC)
  if (ratio == 0)
    if (length(ratio_ind) > accuracy)
      scatter(D_VEC[ratio_ind], DELTA_VEC[ratio_ind], s=size_sc, color="lavender", label=L"$S$%$counter_field: Q-P")
      global counter_field+=1;
    end
  elseif (ratio > 0)
    Burst_p = intersect(Burst_Sync, ratio_ind)
    Spike_p = intersect(Spike_Sync, ratio_ind)
    if (length(Spike_p) > accuracy) || (length(Burst_p) > accuracy)
      if (Spike_p == Burst_p) 
        spike_color = cmap_sp((length(Spike_p) / size)) # Normalize the length to a value between 0 and 1

        scatter(D_VEC[Spike_p], DELTA_VEC[Spike_p], s=size_sc, color=spike_color,label=L"$S$%$counter_field: SS 1:%$ratio")
        global counter_field+=1
      else
        if ((length(Burst_p) > accuracy))
          burst_color = cmap_brst(0.2+(0.5/ratio)) # Normalize the length to a value between 0 and 1

          scatter(D_VEC[Burst_p], DELTA_VEC[Burst_p], s=size_sc, color=burst_color, label=L"$S$%$counter_field: BS 1:%$ratio") 
          global counter_field+=1
        end
        if ((length(Spike_p) > accuracy))
          spike_color = cmap_sp(1.0-(1.0/ratio)) # Normalize the length to a value between 0 and 1

          scatter(D_VEC[Spike_p], DELTA_VEC[Spike_p], s=size_sc, color=spike_color,label=L"$S$%$counter_field: SS 1:%$ratio")
          global counter_field+=1
        end
      end
    end


    Burst_p = intersect(Burst_Sync_ONE, ratio_ind)
    Spike_p = intersect(Spike_Sync_ONE, ratio_ind)
    if (length(Spike_p) > accuracy) || (length(Burst_p) > accuracy)
      if (Spike_p == Burst_p) 
        spike_color = cmap_sp((length(Spike_p) / size)) # Normalize the length to a value between 0 and 1

        scatter(D_VEC[Spike_p], DELTA_VEC[Spike_p], s=size_sc, color=spike_color,label=L"$S$%$counter_field: SS S-P 1:%$ratio")
        global counter_field+=1
      else
        if ((length(Burst_p) > accuracy))
          burst_color = cmap_brst(0.2+(0.5/ratio)) # Normalize the length to a value between 0 and 1

          scatter(D_VEC[Burst_p], DELTA_VEC[Burst_p], s=size_sc, color=burst_color, label=L"$S$%$counter_field: BS S-P 1:%$ratio") 
          global counter_field+=1
        end
        if ((length(Spike_p) > accuracy))
          spike_color = cmap_sp(1.0-(1.0/ratio)) # Normalize the length to a value between 0 and 1

          scatter(D_VEC[Spike_p], DELTA_VEC[Spike_p], s=size_sc, color=spike_color,label=L"$S$%$counter_field: SS S-P 1:%$ratio")
          global counter_field+=1
        end
      end
    end

  elseif (ratio == -1)
    if (length(ratio_ind) > accuracy)
      scatter(D_VEC[ratio_ind], DELTA_VEC[ratio_ind], s=size_sc, color=cmap_dth(0.4), label=L"$S$%$counter_field: D $\varphi_1$")
      global counter_field += 1;
    end
  elseif (ratio == -2)
    if (length(ratio_ind) > accuracy)
      scatter(D_VEC[ratio_ind], DELTA_VEC[ratio_ind], s=size_sc, color=cmap_dth(0.6), label=L"$S$%$counter_field: D $\varphi_2$")
      global counter_field +=1;
    end
  elseif (ratio == -3)
    if (length(ratio_ind) > accuracy)
      scatter(D_VEC[ratio_ind], DELTA_VEC[ratio_ind], s=size_sc, color=cmap_dth(0.8), label=L"$S$%$counter_field: D $\varphi_1, \varphi_2$")
      global counter_field +=1;
    end
  end
end
legend(loc="lower right", fontsize=16, framealpha=1)
# legend(bbox_to_anchor=(1, 1.015), locolor="upper left", fontsize=12)

# xlim(0,  0.4)
title(L"$n_1$ = %$N1, $n_2$ = %$N2,  $\alpha$ = %$alpha_txt", fontsize=20)
xlabel(L"d", fontsize=20)
ylabel(L"\Delta", fontsize=20)