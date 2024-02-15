delta_phi_error = 0.05;

global n1;
global n2;
n1 = 1;
n2 = 1;

global  g1;
global  g2;

global delta;
delta = 0.15;

g1 = 1.01;
g2 = 1.01 + delta;

alpha = 0;

global alpha_text;
alpha_text = '5π/3'; % set alpha value as text

D_CR = [];
global d;
for d = 0:0.001:5
%  for d = 5
  disp(d);
  y0 = [0; 0];
  a=500;
  b=1000;
  opts = odeset('RelTol',2e-13,'AbsTol',1e-100);
  [T,Y] = ode45(@(t,y)eqn(t,y,d, alpha), [a, b], y0, opts);
  [T,Y] = delete_unstb(T,Y, 5);
%   draw(T,Y);
  [OUT, DIFF, ratio] = ms_phase_sync(Y,T);
%   error = 10;
%   OUT = OUT(error:end, 1);
if (ratio == 1)
  if (el_delta(OUT, delta_phi_error))
      D_CR(end+1, 1) = d;
      D_CR(end, 2) = DIFF; 
  end
end
end


plot(D_CR(:,2), D_CR(:,1));




function draw(T,Y)
global g1;
global g2;
global delta;
global d;
global alpha_text;
Y = mod(Y, 2*pi); 
plot(T, Y(:,1),'-b',T, Y(:,2),'-g');
ylim([0 2*pi]); yticks([0 pi 2*pi]); yticklabels(["0" "\pi" "2\pi"]);
xlim([T(1,1) T(end, 1)])

% % Set to show only int

xticks(linspace(T(1,1),T(end, 1), 10));
ax = gca;
ax.XTick = unique(round(ax.XTick));

yticks([0 pi 2*pi]);
xlabel("t");
% xlim([556 T(end,1)])
ylabel('\phi', 'Interpreter','tex');
title(sprintf('d = %g, Δ = %g, α = %s', d, delta, alpha_text));



legend(sprintf('\\gamma_1 = %.2f', g1), ...
       sprintf('\\gamma_2 = %.2f', g2), ...
       'Interpreter', 'tex');

end











function [T,Y] = delete_unstb(T,Y, error)
% error - count of bad spike
% default  error = 5;

  YY = mod(Y, 2*pi);
  spike_start = find(abs(YY(:,2) - 0) < error);
  spike_start = Find_Near_Points(spike_start);
  first_spike = spike_start(error, 1);
  Y = Y(first_spike:end,:);
  T = T(first_spike:end,:);
end


function out = el_delta(OUT, error)
  mn = mean(abs(OUT));
  out =  abs(abs(OUT) - mn) < error;
end



function [OUT, DIFF, ratio] = ms_phase_sync(Y,T)
    YY = mod(Y, 2*pi);

%     DEBUG GRAPHIC


    spike = 2*pi;
    error = 0.05;

%     Number of spikes
    A1 = find(abs(YY(:,1) - spike) < error);
    A2 = find(abs(YY(:,2) - spike) < error);
    A1 = Find_Near_Points(A1);
    A2 = Find_Near_Points(A2);
    if(length(A1) < 1) 
      DIFF = NaN;
      OUT = NaN;
      ratio = 0;
      return;
    end
    
    RATIO = zeros(length(A1(:,1))-2,1);
    for i = 2:length(A1)-1
      RATIO(i-1,1) = sum(A2 < A1(i+1,1)) - sum(A2 < A1(i,1));
    end
    ratio = max(RATIO);

    DIFF = abs(id(A1(1,1),T) - id(A2(1,1),T));


%     FIND NEAR SPIKES
    NEAR = [];

    for i = 1:length(A1)
      [~, index] = max(A2(A2 <= A1(i)));
      if ~isempty(index)
        NEAR(end+1) = A2(index);
      end
    end


    AA1 = id(A1,T);
    NEAR_ID = id(NEAR,T);
    len = length(NEAR_ID(:,1));
    OUT = zeros(len,1);
    AA1 = AA1((end - len + 1):end, 1);
    OUT(:,1) = AA1 - NEAR_ID; 
end


function dy_dt = eqn(~,y,d, alpha)
  global  n1;
  global  n2;
  global  g1;
  global  g2;
  g = [g1; g2];
  no = [n1; n2];
  f = g-sin(y./no);
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
