  using DifferentialEquations
  using LaTeXStrings
  using JLD
  using Statistics
  using Dates
  using Plots
  gr()



  text = L"$\varphi$ нейрон"

  # DELTA = 0.005
  DELTA = 0.01
  # DELTA = 0.05



  N1 = 1;
  N2 = 1;
  name = "fr_dicr_2DIMsin110.01__20240602_1842"

  
  @load "$name.jld2" DATA

  size = length(DATA)
  DATA = reduce(vcat, transpose.(DATA))
  D_VEC = DATA[:,1]
  SIGMA_VEC = DATA[:,2]
  RATIO_S = DATA[:,3]
  RATIO_B = DATA[:,4]

  GOOD = findall(x -> x > 0, RATIO_B[:,1])
  DEAD1 = findall(x -> x == -1.0, RATIO_B[:,1])
  # DEAD2 = findall(x -> x == -2.0, RATIO_B[:,1])
  DEAD3 = findall(x -> x == 0.0, RATIO_B[:,1])

  if length(GOOD) > 0
      ratio_b_good = RATIO_B[GOOD, 1]
      scatter(D_VEC[GOOD], SIGMA_VEC[GOOD], label=L"\Omega_{b}^{1}/\Omega_{b}^{2}", marker_z=ratio_b_good, legend=true)
  end

  if length(DEAD1) > 0
      scatter!(D_VEC[DEAD1], SIGMA_VEC[DEAD1], opacity=0.5, color=colorant"grey44", label=L"D $\; \varphi_1$")
  end
  ratio_b_good = RATIO_B[DEAD3, 1]

  if length(DEAD3) > 0
      ratio_b_good = RATIO_B[DEAD3, 1]
      scatter!(D_VEC[DEAD3], SIGMA_VEC[DEAD3], marker_z=ratio_b_good, legend=true, label=L"D $\varphi_1, \varphi_2$")
  end

  title!(L"$n_1$ = %$N1, $n_2$ = %$N2, Δ = %$DELTA, %$text")
  xlabel!(L"d", guidefontsize=16)
  ylabel!(L"$\sigma$", guidefontsize=16)
  plot!(legendfontsize=10, legend=:topleft)
  # colorbar!(p, label="RATIO_B")
  savefig("$name.png")
