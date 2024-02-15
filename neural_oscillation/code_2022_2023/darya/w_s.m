
global epsilon;
epsilon = 0.002;

d = 0:0.01:1;


n1 = 1;
n2 = 1;

global n
n = [n1; n2];



[w1, w2] = arrayfun(@(x)f1(x),d);

hold on
plot(d, w1,'b');
plot(d, w2, 'g');

xline(d_cr,'-', ' d_{cr}^{s}', 'Interpreter','tex', 'LineWidth', 2);
hold off

legend('w_{s}^{1}','w_{s}^{2}', 'Interpreter','tex');
title('n_1 = n_2 = 1', 'Interpreter','tex');
xlabel("d");
ylabel('w_s', 'Interpreter','tex');

ylim([1 1.4]);

function [w1,w2] = f1(d)
  global n;
  sigma = 1;
  F1_0 = 0;
  F2_0 = 0;
  F = [];
  F(end+1, 1:2) = [0, 0];
  y0 = [pi / 2; pi / 2]; % временные рамки для построения
  a = 2000;
  b = 4000;
  opts = odeset('RelTol', 2e-13, 'AbsTol', 1e-100);
  [T, Y] = ode45(@(t, y)eqn(t, y, F, d, n, sigma), [a, b], y0, opts);

  [A1, A2] = find_spikes(Y);

  w2 = vpa(((length(A2)-1) * 2*pi)/T(end));
  w1 = vpa(((length(A1)-1) * 2*pi)/T(end));
end


function dy_dt = eqn(~, y, F, d, no, sigma)
    yy = mod(y, 2*pi);
    g = [1.01; 1.02];
    f = g - sin(y ./ no);
    exch = d * F(end);
    dy_dt = f - exch;
    Fo = chech_condition(yy, sigma);
    F(end+1,1:2) = [Fo(1), Fo(2)];
    disp(Fo);
end



function F = chech_condition(y, sigma)
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

  F = [F2, F1];
end

function [A1, A2] = find_spikes(Y)
  error = 0.05;
  spike = 2 * pi;
  YY = mod(Y, 2 * pi);
  A1 = find(abs(YY(:, 1) - spike) < error);
  A2 = find(abs(YY(:, 2) - spike) < error);
  A1 = Find_Near_Points(A1);
  A2 = Find_Near_Points(A2);
end

function AA1 = Find_Near_Points(A1)
  AA1 = A1;
  tolerance = 2;
  diffs = diff(A1);
  indices_to_remove = find(diffs < tolerance);
  AA1(indices_to_remove + 1) = [];
end


function tilda = tld(n)
    tilda = n*pi/2 - floor(n/4)*2*pi;
    if abs(tilda - pi/2) < 0.02
        tilda = 3*pi/2;
    else
        tilda = abs(tilda - pi);
    end
end

function AA1 = id(AA1, T)
    AA1 = num2cell(AA1);
    AA1 = T([AA1{:}], :);
end