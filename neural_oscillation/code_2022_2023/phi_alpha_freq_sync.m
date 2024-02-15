% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  Freq synchronization modes for phi-alpha system for 2-neurons
%
%  alpha - offset parameter of system
%  d - sync parameter of system
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
addpath('additional\progress_bar\'); % add progress bar

format longG
n = [3; 3];

delta_max = 0.01;

g1 = 1.01;
g = [g1; g1 + delta_max];

alpha = 2*pi/3;
color = 'green';
global alpha_text
alpha_text = '2π/3';

DiagOut = diag(alpha, g, n);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function  DiagOut = diag(alpha, g, n)
  d_list = 0:0.0001:0.09;
  DiagOut = zeros(length(d_list), 3);
  parfor m = 1:length(d_list)
    d = d_list(m);
    y0 = [0; 0];
    a=2000;
    b=20000;
    opts = odeset('RelTol',2e-13,'AbsTol',1e-100, 'Refine', 10);
    [T,Y] = ode45(@(t,y)eqn(g, n, t, y, d, alpha), [a, b], y0, opts);
% draw(T,Y,g,d);
    [freq_ratio, err] = FREQ(T,Y, n);
    DiagOut(m, :) = [freq_ratio, d, err];
  end
end

function [ratio, err] = FREQ(T,Y, n)

  TT = ms(Y,T,n);
  TT1 = nonzeros(TT(:,1));
  TT2 = nonzeros(TT(:,2));

  if ((length(TT1) < 1) && (length(TT2) < 1))
      err = 3;
  elseif (length(TT2) < 1)
      err = 2;
  elseif (length(TT1) < 1)
      err = 1;
  else
      err = 0;
  end

  if (err > 0)
      ratio = 0;
      return;
  end

% FOR BURST FREQ USE THIS BLOCK
% TT1 = sum(TT1);
% TT2 = sum(TT2);
  TT1 = del_middle_burst_int(n(1), TT1);
  TT2 = del_middle_burst_int(n(2), TT2);
  TT1 = mean(TT1);
  TT2 = mean(TT2);

  w1 = vpa(2*pi/TT1);
  w2 = vpa(2*pi/TT2);
  ratio = w1/w2;
end

function T_DIFF = del_middle_burst_int(n, T)
  if (n > 1)
    T = sort(T,'descend');
    T_DIFF = T(2:size(T));
  else 
    T_DIFF = T;
  end
end

function OUT = ms(Y,T,n)
  OUT = zeros(max(n),length(n));
  for ind = 1:length(n)
      UNSTBL = delete_unstb(T,Y, n, 20, ind);
      for k = 1:n(ind)
          [A,err] = find_spikes(Y, n, ind);
          if ((err) > 0)
              return;
          end
          start = find(abs((A - UNSTBL)) < 10);
          if (isnan(UNSTBL))
              start = 1;
          end
          Diff_points = A(start:end, 1);
          p1 = Diff_points(k);
          p2 = Diff_points(k+1);
          OUT(k, ind) = id(p2,T) - id(p1,T);
      end
  end
end

function AA1 = Find_Near_Points(A1)
  AA1 = A1;
  tolerance = 2;
  diffs = diff(A1);
  indices_to_remove = find(diffs < tolerance);
  AA1(indices_to_remove + 1) = [];
end

function AA1 = id(AA1, T)
    AA1 = num2cell(AA1);
    AA1 = T([AA1{:}], :);
end

function UNSTBL = delete_unstb(T, Y, n, error, ind)
    UNSTBL = NaN;
    [A, err] = find_spikes(Y, n, ind);
    if ((error == 0) && (err > 0))
        return;
    end
    num = n(ind);
    B = zeros(floor(size(A,1)/num),1);
    for m = 1:num:length(A)
        B(ceil(m/num),1) = A(m,1);
    end
    if ((length(B)) < (error + 2))
        return;
    end
    UNSTBL = B(error, 1);
end

function [A, err] = find_spikes(Y, n, ind)
  tilda = tld(n(ind));
  error = 0.05;
  YY = mod(Y, 2 * pi);
  A = find(abs(YY(:, ind) - tilda) < error);
  A = Find_Near_Points(A);

  if (find(abs(YY(:, ind) - 2*pi) < error))
    err = 0;
  else
    err = 1;
  end
end

function tilda = tld(n)
  tilda = n*pi/2 - floor(n/4)*2*pi;
  if abs(tilda - pi/2) < 0.02
      tilda = 3*pi/2;
  else
      tilda = abs(tilda - pi);
  end
end

function dy_dt = eqn(g, n, ~, y, d, alpha)
    f = g - sin(y ./ n);
    exch = [d * sin(y(2) - y(1) - alpha); d * (sin(y(1) - y(2) - alpha))];
    dy_dt = f + exch;
end

function draw(T, Y, g, d)
    global alpha_text;
  delta_max = g(2) - g(1);
  Y = mod(Y, 2 * pi);
  plot(T, Y(:, 1), '-b', T, Y(:, 2), '-g');
  ylim([0 2 * pi]); yticks([0 pi 2 * pi]); yticklabels(["0" "\pi" "2\pi"]);
  ax = gca;
  ax.XTick = unique(round(ax.XTick));
  yticks([0 pi 2 * pi]);
  xlabel("t");
  % xlim([556 T(end,1)])
  ylabel('\phi', 'Interpreter', 'tex');
  title(sprintf('d = %g, Δ = %g, α = %s', d, delta_max, alpha_text));
  legend(sprintf('\\gamma_1 = %.2f', g(1)), ...
      sprintf('\\gamma_2 = %.2f', g(2)), ...
      'Interpreter', 'tex');
end