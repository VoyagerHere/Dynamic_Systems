using Plots

d =  0.0821
sigma = 1.58
n1 = 1
n2 = 1
n = [n1, n2]
a = 0
b = 1000
h = 1/100
T = a:h:(b-h)
Y = zeros(length(T), 2)
y0 = [pi/2, pi/2]
Y[1,:] = y0
global F = [0, 0]

function eqn(t, y, d, no, sigma)
    global F
    yy = mod.(y, 2*pi)
    g = [1.001, 1.002]
    f = g - sin.(y ./ no)
    exch = d * [F[2], F[1]]
    dy_dt = f - exch
    F1, F2 = chech_condition(yy, sigma)
    F = [F1, F2]
    return dy_dt
end

function chech_condition(y, sigma)
    F1 = (y[1] > pi/2 - sigma) && (y[1] < pi/2 + sigma) ? 0 : 1
    F2 = (y[2] > pi/2 - sigma) && (y[2] < pi/2 + sigma) ? 0 : 1
    return F1, F2
end

for i = 1:(length(T)-1)
    dy_dt = eqn(T[i], Y[i,:], d, n, sigma)
    Y[i+1,1] = Y[i,1] + dy_dt[1]*(h)
    Y[i+1,2] = Y[i,2] + dy_dt[2]*(h)
end

Y = mod.(Y, 2*pi)

function draw(T, Y, d, sigma)
    p = plot(T, Y[:, 1], line = (:red), label = "\$\\phi_1\$")
    plot!(T, Y[:, 2], line = (:darkgreen), label = "\$\\phi_2\$")
    ylims!(0, 2 * pi)
    yticks!([0, pi, 2 * pi], ["0", "\\pi", "2\\pi"])
    xlabel!("t")
    ylabel!("\$\\phi\$")
    title!("d = $d, σ = $sigma")
    display(p)
end

draw(T, Y, d, sigma)