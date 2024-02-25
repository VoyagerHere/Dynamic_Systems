using DifferentialEquations
using Plots
using LaTeXStrings
using JLD
using Statistics
using Dates


const NUM_OF_COMPUTE_RES = 4;
G1 = 1.01;
G2 = 1.02;
alpha_txt = "π"



@load "pi_8_2_2.jld2" DATA SYNC

size = length(DATA)
DATA = reduce(vcat,transpose.(DATA))
D = DATA[:,1]
RATIO = round.(DATA[:,2], digits=1, RoundUp)
DELTA = DATA[:,3]

unique_ratios = unique(RATIO)
for ratio in unique_ratios
  if (ratio == NaN)
    continue
  end
    

end

SYNC_BS = SYNC[:, 1]
SYNC_SP = SYNC[:, 2]

plot(D, RATIO,)
title!(L"\alpha = %$alpha_txt, \gamma_{1}=%$G1, \gamma_{2}=%$G2")
ylims!(0,  1)
xlabel!(L"d")
ylabel!(L"\frac{w_1}{w_2}")