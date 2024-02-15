n1 = 1;
n2 = 1;
n = [n1; n2];
delta_phi_error = 0.05;
g1 = 1.01;

global delta_max
delta_max = 0.001;
global alpha_text
alpha_text = '2π/3';
alpha = 2*pi/3;

% Table with gamma
g_list = linspace(1.01, g1 + delta_max, 1);
DIFF = [];
DIFF = gamma_iterate(g1, n, g_list, DIFF, alpha);

plot(DIFF);


function DIFF = gamma_iterate(g1, n, g_list, DIFF, alpha)
  delta_phi_error = 0.005;
  for k = 1 : length(g_list)
    g2 = g_list(k);
    DIFF = d_iterate([g1; g2], n, alpha, delta_phi_error);
  end
end

% function DIFF = gamma_iterate(g1, n, g_list, GAMMA_TABLE, alpha)
%   delta_phi_error = 0.005;
%   for k = 1 : length(g_list)
%     g2 = g_list(k);
%     DIFF = d_iterate([g1; g2], n, alpha, delta_phi_error);
%   end
% end


function  DIFF = d_iterate(g, n, alpha, delta_phi_error)
    d = 0.012;  % for 2pi/3
    disp(g(2));
    y0 = [0; 0];
    a=1500;
    b=10000;
%     b=2500; % for pi/8
    if((g(1) == g(2))&&(d > 0.026))
      b = 100000;
    end
    opts = odeset('RelTol',2e-13,'AbsTol',1e-100, 'Refine', 10);
    [T,Y] = ode45(@(t,y)eqn(g, n, t,y,d, alpha), [a, b], y0, opts);
    [A1, A2, err] = find_spikes(Y);
    if (err == 0)
      [T,Y] = delete_unstb(T,Y, A1, A2, 10);
    end
    [DIFF, diff, ratio] = ms_phase_sync(T,Y);
    if (el_delta(abs(DIFF), delta_phi_error)) % phase sync
      disp('YEAS');
%     draw(T,Y, DIFF, g, d); % DEBUG
    end
end

function [A1, A2, err] = find_spikes(Y)
% err = 1  - Death
  error = 0.05;
  spike = 2*pi;
  YY = mod(Y, 2*pi);
  A1 = find(abs(YY(:,1) - spike) < error);
  A2 = find(abs(YY(:,2) - spike) < error);
  A1 = Find_Near_Points(A1);
  A2 = Find_Near_Points(A2);
  if(length(A1) < 1)
    if(length(A2) < 1)
      err = 2;
    else
      err = 1;
    end
  else 
    err = 0;
  end
end

function [T,Y] = delete_unstb(T,Y, A1, A2, error)
% error - count of bad spike
% default  error = 5;
%   draw(T,Y,0);
  error = min(error, max((length(A1)-4),0));
  if (error == 0)
      return;
  end

  first_spike = A1(error, 1);
  Y = Y(first_spike:end,:);
  T = T(first_spike:end,:);
end

function out = el_delta(DIFF, error)
  mn = mean(abs(DIFF));
  out = abs(abs(DIFF) - mn) < error;
end

function [DIFF, diff, ratio] = ms_phase_sync(T,Y)
    [A1, A2, err] = find_spikes(Y);
    if (err > 0) % death
       DIFF = NaN;
       diff = NaN;
       ratio = -err;
       return;
    end

    RATIO = zeros(length(A1(:,1))-2,1);
    for i = 2:length(A1)-1
      RATIO(i-1,1) = sum(A2 < A1(i+1,1)) - sum(A2 < A1(i,1));
    end

    ratio = mode(RATIO);    
    NEAR = [];  % find near spikes
    for i = 1:length(A1)
      [~, index] = max(A2(A2 <= A1(i)));
      if ~isempty(index)
        NEAR(end+1) = A2(index);
      end
    end
    AA1 = id(A1,T);
    NEAR_ID = id(NEAR,T);
    len = length(NEAR_ID(:,1));
    DIFF = zeros(len,1);
    AA1 = AA1((end - len + 1):end, 1);
    DIFF(:,1) = AA1 - NEAR_ID;
    diff = DIFF(end, 1);
end

function draw(T,Y, DIFF, g, d)
    global delta_max;
    global alpha_text;
  Y = mod(Y, 2*pi); 
  plot(T, Y(:,1),'-b',T, Y(:,2),'-g');
  ylim([0 2*pi]); yticks([0 pi 2*pi]); yticklabels(["0" "\pi" "2\pi"]);
  xlim([T(1,1) T(end, 1)])
  xticks(linspace(T(1,1),T(end, 1), 10)); % set to show only int
  ax = gca;
  ax.XTick = unique(round(ax.XTick));
  yticks([0 pi 2*pi]);
  xlabel("t");
  % xlim([556 T(end,1)])
  ylabel('\phi', 'Interpreter','tex');
  title(sprintf('d = %g, Δ = %g, α = %s, |φ₁ - φ₂| ≈ %.2f', d, delta_max, alpha_text, DIFF));
  legend(sprintf('\\gamma_1 = %.2f', g(1)), ...
          sprintf('\\gamma_2 = %.2f', g(2)), ...
          'Interpreter', 'tex');
end

function dy_dt = eqn(g, n, ~, y,d, alpha)
  f = g-sin(y./n);
  exch = [d * sin(y(2) - y(1) - alpha); d * (sin(y(1) - y(2) - alpha))];
  dy_dt = f+exch;
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
