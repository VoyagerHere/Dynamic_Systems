using DifferentialEquations
using LaTeXStrings
using JLD
using Statistics
using Dates

using PythonPlot
pygui(true)

text = L"$\varphi$ нейрон"

# DELTA = 0.005
# DELTA = 0.01
DELTA = 0.05

N1 = 3
N2 = 3

name = "fr_dicr_2DIMsin330.05__20240609_0623.jld2"

@load "$name" DATA

PythonPlot.matplotlib.rcParams["font.size"] = 14
size_sc = 4

size = length(DATA)
DATA = reduce(vcat, transpose.(DATA))
D_VEC = DATA[:,1]
SIGMA_VEC = DATA[:,2]
RATIO_S = DATA[:,3]
RATIO_B = DATA[:,4]

GOOD = findall(x -> x > 0, RATIO_B[:,1])
DEAD1 = findall(x -> x == -1.0, RATIO_B[:,1])
DEAD3 = findall(x -> x == 0.0, RATIO_B[:,1])

# Create the figure and axis objects
fig, ax = subplots()

if length(GOOD) > 0
    ratio_b_good = RATIO_B[GOOD, 1]
    # Normalize the ratio_b_good values to [0, 1]
    # normalized_ratio_b_good = (ratio_b_good .- minimum(ratio_b_good)) ./ (maximum(ratio_b_good) - minimum(ratio_b_good))
    scatter_plot = ax.scatter(D_VEC[GOOD], SIGMA_VEC[GOOD], s=size_sc, c=ratio_b_good, cmap="viridis", label=L"\Omega_{b}^{1}/\Omega_{b}^{2}")
    cbar = fig.colorbar(scatter_plot, ax=ax)
    cbar.set_label(L"$\Omega_{b}^{1}/\Omega_{b}^{2}$")
end

# if length(DEAD1) > 0
#     ax.scatter(D_VEC[DEAD1], SIGMA_VEC[DEAD1], s=size_sc, alpha=0.5, color="grey", label=L"D $\; \varphi_1$")
# end

# if length(DEAD3) > 0
#     ax.scatter(D_VEC[DEAD3], SIGMA_VEC[DEAD3], s=size_sc, label=L"D $\varphi_1, \varphi_2$")
# end

# ax.legend(loc="lower right", fontsize=16, framealpha=1)

ax.set_title(L"$n_1$ = %$N1, $n_2$ = %$N2, Δ = %$DELTA, %$text")
ax.set_xlabel(L"d", fontsize=16)
ax.set_ylabel(L"\sigma", fontsize=20)