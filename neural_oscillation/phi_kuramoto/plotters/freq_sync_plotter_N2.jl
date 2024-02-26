using DifferentialEquations
using LaTeXStrings
using JLD
using Statistics
using Dates
import PythonPlot as pyplot
pyplot.pygui(true)

const NUM_OF_COMPUTE_RES = 4;
G1 = 1.01;
G2 = 1.02;
name = ""
alpha_txt = "π/8"
N1 = 2;
N2 = 2;

@load "untitled__20240226_1635.jld2" DATA

size = length(DATA)
DATA = reduce(vcat,transpose.(DATA))
D_VEC = DATA[:,1]
W_1 = DATA[:,2]
W_2 = DATA[:,3]
RATIO_W = DATA[:,4]


RATIO_W = round.(RATIO_W; digits = 3);

GOOD =  findall(x -> x > 0, RATIO_W[:,1])
DEAD1 = findall(x -> x ==  -1.0, RATIO_W[:,1])
DEAD2 = findall(x -> x ==  -2.0, RATIO_W[:,1])
DEAD2 = findall(x -> x ==  -3.0, RATIO_W[:,1])

pyplot.plot(D_VEC[GOOD], RATIO_W[GOOD])

x = range(0,  1, length(D_VEC))
y = sin.(x)
y = range(0,  1, length(D_VEC))

pyplot.fill_between(D_VEC[DEAD1], x[DEAD1], color="blue", alpha=0.5)
pyplot.show()
#plot(D_VEC[DEAD2], x)
#plot(D_VEC[DEAD3], x)

title(L"$n_1$ = %$N1, $n_2$ = %$N2,  $\alpha$ = %$alpha_txt, $\gamma_{1}$=%$G1, $\gamma_{2}$=%$G2")
ylims(0,  1)
xlabel(L"d")
ylabel(L"\frac{w_1}{w_2}")