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
DEATH_MATR = Int.(reduce(vcat,transpose.(DEATH)))

D_VEC_1 = DATA[:,1]
D_VEC_2 = DATA[:,2]
RATIO_VEC_12 = DATA[:,3]
RATIO_VEC_23 = DATA[:,4]
DELTA_VEC = DATA[:,5]

unique_ratios = unique(RATIO_VEC)


counter_field = 1;

PythonPlot.matplotlib.rcParams["font.size"] = 14
size_sc = 4

DEAD1 = findall(x -> x ==  1, DEAD[:,1])
DEAD2 = findall(x -> x ==  1, DEAD[:,2])
DEAD3 = findall(x -> x ==  1, DEAD[:,3])
DEAD_ALL = intersect(DEAD1, DEAD2, DEAD3)
QUASI_PERIODIC = findall(x -> x ==  0, DEAD_ALL)

global counter_field+=1;

if(length(DEAQUASI_PERIODICD1) > accuracy)
  scatter(D_VEC_1[QUASI_PERIODIC], DELTA_VEC[QUASI_PERIODIC], color = "snow", s=size_sc, label=L"$S_{1}$: Q-P")
  global counter_field+=1;
end

if(length(DEAD1) > accuracy)
  scatter(D_VEC_1[DEAD1], DELTA_VEC[DEAD1], s=size_sc, color = "slategrey", label=L"$S_{counter_field}$: D $\varphi_1$")
  global counter_field+=1;
end

if(length(DEAD2) > accuracy)
  scatter(D_VEC_1[DEAD2], DELTA_VEC[DEAD2], s=size_sc, color = "darkgrey", label=L"$S_{1}$: D $\varphi_2$")
  global counter_field+=1;
end

if(length(DEAD3) > accuracy)
  scatter(D_VEC_1[DEAD3], DELTA_VEC[DEAD3], s=size_sc, color = "dimgrey", label=L"$S_{1}$: D $\varphi_3$")
  global counter_field+=1;
end

if(length(DEAD_ALL) > accuracy)
  scatter(D_VEC_1[DEAD_ALL], DELTA_VEC[DEAD_ALL], s=size_sc, color = "black", label=L"$S_{1}$: D $\varphi_1, \varphi_2, \varphi_3$")
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