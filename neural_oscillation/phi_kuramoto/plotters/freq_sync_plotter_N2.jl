using DifferentialEquations
using LaTeXStrings
using JLD
using Statistics
using Dates
using Plots
gr()

alpha_txt = "π/2"
N1 = 3;
N2 = 3;
@load "pi_1_2__20240403_1505.jld2" DATA

size = length(DATA)
DATA = reduce(vcat,transpose.(DATA))
D_VEC = DATA[:,1]
RATIO_S = DATA[:,2]
RATIO_B = DATA[:,3]


GOOD =  findall(x -> x > 0, RATIO_B[:,1])
DEAD1 = findall(x -> x ==  -1.0, RATIO_B[:,1])
DEAD2 = findall(x -> x ==  -2.0, RATIO_B[:,1])
DEAD3 = findall(x -> x ==  -3.0, RATIO_B[:,1])


rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])


if length(GOOD) > 0
  plot(D_VEC[GOOD], RATIO_B[GOOD], label=L"\Omega_{b}^{1}/\Omega_{b}^{2}", c = "aqua")
  scatter!([minimum(D_VEC)],[0], label=" ", ms=0, mc=:white, msc=:white)
end
if length(DEAD1) > 0
  plot!(rectangle(maximum(D_VEC[DEAD1])-minimum(D_VEC[DEAD1]),1,minimum(D_VEC[DEAD1]),0), opacity=.5, color = colorant"grey44", label = L"D $\; \varphi_1$")
  scatter!([minimum(D_VEC)],[0], label=" ", ms=0, mc=:white, msc=:white)
end
if length(DEAD2) > 0
  plot!(rectangle(maximum(D_VEC[DEAD2])-minimum(D_VEC[DEAD2]),1,minimum(D_VEC[DEAD2]),0), opacity=.5, color = colorant"grey34", label = L"D $\varphi_2$")
  scatter!([minimum(D_VEC)],[0], label=" ", ms=0, mc=:white, msc=:white)
end
if length(DEAD3) > 0
  plot!(rectangle(maximum(D_VEC[DEAD3])-minimum(D_VEC[DEAD3]),1,minimum(D_VEC[DEAD3]),0), opacity=1, color = colorant"black",  label = L"D $\varphi_1, \varphi_2$")
end

# # Calculate center points for each area
# center_DEAD1 = (minimum(D_VEC[DEAD1]) + maximum(D_VEC[DEAD1])) /  2
# center_DEAD2 = (minimum(D_VEC[DEAD2]) + maximum(D_VEC[DEAD2])) /  2
# center_DEAD3 = (minimum(D_VEC[DEAD3]) + maximum(D_VEC[DEAD3])) /  2

# annotate!(center_DEAD1,  0.5, text("DEAD1", :center,  5))
# annotate!(center_DEAD3,  0.5, text("Death3", :center,  5))

title!(L"$n_1$ = %$N1, $n_2$ = %$N2, $\alpha$ = %$alpha_txt")
ylims!(0.7,  1)
xlabel!(L"d", guidefontsize=16)
ylabel!(L"\Omega_{1}/\Omega_{2}", guidefontsize=16)
plot!(legendfontsize=10, legend=:outertopright) 
savefig("freq_plot.png")