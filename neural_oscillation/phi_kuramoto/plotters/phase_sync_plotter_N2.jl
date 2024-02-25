using DifferentialEquations
using Plots
using LaTeXStrings
using JLD
using Statistics
using Dates


const NUM_OF_COMPUTE_RES = 4;
G1 = 1.01;
G2 = 1.02;
name = ""
alpha_txt = "π"

@load "$name.jld2" DATA, SYNC

size = length(DATA)
D = DATA[:,1]
RATIO = DATA[:,2]
DELTA = DATA[:,3]

SYNC_BS = SYNC[:, 1]
SYNC_SP = SYNC[:, 2]

plot(D, RATIO,)
title!(L"\alpha = %$alpha_txt, \gamma_{1}=%$G1, \gamma_{2}=%$G2")
ylims!(0,  1)
xlabel!(L"d")
ylabel!(L"\frac{w_1}{w_2}")