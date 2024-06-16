d =  0.001;
sigma = 0;

n1 = 1;
n2 = 1;
n = [n1; n2];

a = 0;
b = 6000;
h = 1/100;
T = a:h:(b-h);
T = transpose(T);


Y(length(T), 2) = 0;
y0 = [pi/2; pi/2];
Y(1,1:2) = y0;
F1_0 = 0;
F2_0 = 0;
F = [F1_0; F2_0];

% Solve ODE for each time step and update F
for i = 1:((length(T)-1))
    [dy_dt, F] = eqn(T(i), Y(i,:), d, n, sigma, F);
    Y(i+1,1) = Y(i,1) + dy_dt(1)*(h);
    Y(i+1,2) = Y(i,2) + dy_dt(2)*(h);
end
Y = mod(Y, 2*pi);
draw(T, Y, d, sigma);


function [dy_dt, F] = eqn(~, y, d, no, sigma, F)
    y = transpose(y);
    yy = mod(y, 2*pi);
    g = [1.001; 1.002];
    f = g - cos(y ./ no);
    exch = d * [F(2); F(1)];
    dy_dt = f - exch;
    [F1, F2] = chech_condition(yy, sigma);
    F = [F1; F2];
end


function [F1, F2] = chech_condition(y, sigma)
  if ((y(1) > pi/2 - sigma) && (y(1) < pi/2 + sigma)) 
    F1 = 0;
  else
    F1 = 1;
  end

  if ((y(2) > pi/2 - sigma) && (y(2) < pi/2 + sigma)) 
    F2 = 0;
  else
    F2 = 1;
  end
end

function draw(T, Y, d, sigma)
    plot(T, Y(:, 1), '-r', T, Y(:, 2), '-b');
    legend('\phi_1', '\phi_2')
    ylim([0 2 * pi]);
    yticks([0 pi 2 * pi]);
    yticklabels(["0" "\pi" "2\pi"]);
    xlabel("t");
    ylabel('\phi', 'Interpreter', 'tex');
    title(sprintf('d = %g, σ = %g', d, sigma));
end