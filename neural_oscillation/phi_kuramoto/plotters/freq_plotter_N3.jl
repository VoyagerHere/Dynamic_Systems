using DifferentialEquations
using LaTeXStrings
using JLD
using Statistics
using Dates
using Plots
gr()


alpha_txt = "2π/3"
N1 = 3;
N2 = 3;
N3 = 3;
name = "\Omega_{1,2}/\Omega_{2,3}"
@load "pi_2_3__3_3_3__20240327_0602.jld2" DATA W DEATH
 
size = length(DATA)
DATA = reduce(vcat,transpose.(DATA))
DEAD = reduce(vcat,transpose.(DEATH))
W_ARR = reduce(vcat,transpose.(W))

D_VEC_1 = DATA[:,1]
D_VEC_2 = DATA[:,2]


D_VEC = D_VEC_1

W_B_1_2 = W_ARR[:,2]
W_B_1 = W_B_1_2[1]
W_B_2 = W_B_1_2[2]
W_B_2_3 = W_ARR[:,4]
W_B_3 = W_B_2_3[2]



DEAD1 = findall(x -> x ==  1, DEAD[:,1])
DEAD2 = findall(x -> x ==  1, DEAD[:,2])
DEAD3 = findall(x -> x ==  1, DEAD[:,3])


rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])



plot(D_VEC, W_B_1, label=L"\Omega_{b}^{1}", c = "azure")
scatter!([minimum(D_VEC)],[0], label=" ", ms=0, mc=:white, msc=:white)


plot(D_VEC, W_B_2, label=L"\Omega_{b}^{2}", c = "aqua")
scatter!([minimum(D_VEC)],[0], label=" ", ms=0, mc=:white, msc=:white)

plot(D_VEC, W_B_3, label=L"\Omega_{b}^{3}", c = "aqua")
scatter!([minimum(D_VEC)],[0], label=" ", ms=0, mc=:white, msc=:white)




if length(DEAD1) > 0
  plot!(rectangle(maximum(D_VEC[DEAD1])-minimum(D_VEC[DEAD1]),1,minimum(D_VEC[DEAD1]),0), opacity=.5, color = colorant"grey24", label = L"D $\varphi_1$")
  scatter!([minimum(D_VEC)],[0], label=" ", ms=0, mc=:white, msc=:white)
end
if length(DEAD2) > 0
  plot!(rectangle(maximum(D_VEC[DEAD2])-minimum(D_VEC[DEAD2]),1,minimum(D_VEC[DEAD2]),0), opacity=.5, color = colorant"grey34", label = L"D $\varphi_1, \varphi_2$")
  scatter!([minimum(D_VEC)],[0], label=" ", ms=0, mc=:white, msc=:white)
end
if length(DEAD3) > 0
  plot!(rectangle(maximum(D_VEC[DEAD3])-minimum(D_VEC[DEAD3]),1,minimum(D_VEC[DEAD3]),0), opacity=1, color = colorant"black",  label = L"D $\varphi_1, \varphi_2, \varphi_3$")
end

# # Calculate center points for each area
# center_DEAD1 = (minimum(D_VEC[DEAD1]) + maximum(D_VEC[DEAD1])) /  2
# center_DEAD2 = (minimum(D_VEC[DEAD2]) + maximum(D_VEC[DEAD2])) /  2
# center_DEAD3 = (minimum(D_VEC[DEAD3]) + maximum(D_VEC[DEAD3])) /  2

# annotate!(center_DEAD1,  0.5, text("DEAD1", :center,  5))
# annotate!(center_DEAD3,  0.5, text("Death3", :center,  5))


title!(L"$n_1$ = %$N1, $n_2$ = %$N2, $n_3$ = %$N3, $\alpha$ = %$alpha_txt")
# ylims!(0,  1)
#  xlims!(0,  0.03)
xlabel!(L"d", guidefontsize=16)
ylabel!(L"%$name", guidefontsize=16)
plot!(legendfontsize=10, legend=:outertopright) 
savefig("freq_plot.png")