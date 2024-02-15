format longG

% n 
global n1;
global n2;
n1 = 1;
n2 = 1;

% g
global  g1;
global  g2;
g1 = 1.01;

global delta_max;
delta_max = 0.01;

% alpha
color = 'green';
alpha = 2*pi/3;
alpha_text = '2π/3'; % set alpha value as text




DiagOut_1 = diag(alpha);

% Close plot of phi
close;

% Find there ratio is changed
% idx = find(diff(RATIO(:,1)) ~= 0) + 1;
% idx = idx - 1;



% % PLOT
hold on

% % ROUND
DiagOut_1(:,1) = round(DiagOut_1(:,1), 3);

% % Get unique y-values
y_values = unique(DiagOut_1(:,1)); 

% % Crop y_values to the element before the first NaN element
nan_index = find(isnan(y_values), 1, 'first');
if ~isempty(nan_index)
    y_values = y_values(1:nan_index-1);
end



%     Find Death
death_x = Find_Death(DiagOut_1);
area(death_x, ones(size(death_x)));

for i = 1:length(y_values)
    value = y_values(i);
    if (isnan(value)) 
      continue;
    end
    x_values = DiagOut_1(DiagOut_1(:,1) == value, 2); % Get x-values where y = y_values(i)
    plot(x_values, y_values(i)*ones(size(x_values)), 'b', 'LineWidth', 2); % Plot bold and blue horizontal line at y = y_values(i)    hold on;
end


% legend('$\alpha = \frac{\pi}{6}$', 'Interpreter', 'latex');
hold off
xlim([min(DiagOut_1(:, 2)), max(DiagOut_1(:, 2))]);
title(sprintf('Δ = %g, α = %s', delta_max, alpha_text));
h = ylabel('$\omega_1 / \omega_2$', 'Interpreter', 'latex');
set(h, 'VerticalAlignment', 'bottom');
xlabel('d');



function nan_x_values = Find_Death(DiagOut_1)
  y_values = DiagOut_1(:,1); % Get all y-values
  nan_indices = isnan(y_values); % Find indices of NaN elements
  nan_x_values = DiagOut_1(nan_indices, 2); % Get x-values corresponding to NaN y-values
end



function  DiagOut = diag(alpha)
  % delta
  delta_phi_error = 0.005;
  global delta_max;

  global g1;
  global g2;
  gg2 = linspace(g1, g1 + delta_max, 1);
  
  DiagOut = [];
  for k = 1 : length(gg2)
    g2 = gg2(k);
    Delta = g2 - g1;
    disp(Delta);
    for d = 0:0.00001:1.0
      disp(d);
      y0 = [0; 0];
      a=500;
      b=1000;
      opts = odeset('RelTol',2e-13,'AbsTol',1e-100);
      [T,Y] = ode45(@(t,y)eqn(t,y,d, alpha), [a, b], y0, opts);
%       [T,Y] = delete_unstb(T,Y, 5);
      freq_ratio = FREQ(T,Y);

      len = size(DiagOut,1);

      DiagOut(len+1, 1) = freq_ratio;
      DiagOut(len+1, 2) = d;
      
%       GRAPHIC

%       Y = mod(Y, 2*pi); 
%       plot(T, Y(:,1),'-b',T, Y(:,2),'-g');
%       ylim([0 2*pi]); yticks([0 pi 2*pi]); yticklabels(["0" "\pi" "2\pi"]);
    end
  end
end

function ratio = FREQ(T,Y)

%   ! FIX FOR n > 1
%   global n1;
%   global n2;

  [TT2] = ms(Y,T,2);
  [TT1] = ms(Y,T,1);

  TT2 = mean(TT2);
  TT1 = mean(TT1);

  w2 = 2*pi/TT2;
  w1 = 2*pi/TT1;
  ratio = w1/w2;
end

function dy_dt = eqn(~,y,d, alpha)
  global  n1;
  global  n2;
  global  g1;
  global  g2;
  g = [g1; g2];
  no = [n1; n2];
  f = g-sin(y./no);
  exch = [d * sin(y(2)-y(1) - alpha); d * (sin(y(1) - y(2) - alpha))];
%   exch = [d;-d]*sin(y(2)-y(1) - alpha);
  dy_dt = f+exch;
end

function OUT = ms(Y,T,ind)
    YY = mod(Y(:,ind), 2*pi);
    spike = 2*pi;
    error = 0.01;

    A1 = find(abs(YY - spike) < error);
    A1 = Find_Near_Points(A1);
%     spike_count = length(A1);
    AA1 = id(A1,T);
    OUT = diff(AA1);
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