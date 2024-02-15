warning('off');


d = 0:0.001:1.3;
% d = 1:0.1:1.3;

% d = [1, 1.2];
% d = 0;

global epsilon;
epsilon = 0.002;
global D;
D = [];
global n1;
global n2;
n1 = 3;
n2 = 6 ;


[w1, w2] = arrayfun(@(x)f1(x),d);
hold on
plot(d, w1,'g');
plot(d, w2, 'b');

d_cr = find_cr(D);

xline(d_cr,'-', ' d_{cr}^{b}', 'Interpreter','tex', 'LineWidth', 2);
hold off

legend('w_{b}^{1}','w_{b}^{2}', 'Interpreter','tex');
title('n_1 = 4, n_2 = 6', 'Interpreter','tex');
xlabel("d");
ylabel('w_b', 'Interpreter','tex');

function [w1, w2] = f1(d)
  global D;
  global epsilon;
  y0 = [0; 0];
  a=0;
  b=10000;
  opts = odeset('RelTol',2e-13,'AbsTol',1e-100);
  [T,Y] = ode45(@(t,y)eqn(t,y,d), [a, b], y0, opts);
  
  w1 = vpa((Y(end,1) - Y(502,1))/(b*2*pi));
  w2 = vpa((Y(end,2) - Y(502,2))/(b*2*pi));

  if (abs(w2 - w1) < epsilon)
    len = size(D,1);
    D(len+1,1) = d;
    D(len+1,2) = w1;
  end
end

function d_cr = find_cr(D)
    error = 0.001;
    error_int = 0.3;
    w_med = D(1,2);
    d_cr = D(1,1);
    delta = 0;
    for i = 1:length(D(:,1)) 
%       length of interval > error_int
        if (delta > error_int)
            break;
        end
%       find where distance between new w and w prev < error
%       and update interval from d_cr to new d
        if (abs(D(i,2) - w_med) < error)
            delta = D(i,1) - d_cr;
        else
            delta = 0;
            d_cr = D(i,1);
            w_med = D(i,2);
        end
    end
end

function dy_dt = eqn(~,y,d)
  global n1;
  global n2;
  g = [1.01; 1.02];
  no = [n1; n2];
  f = g-sin(y./no);
  exch = [d;-d]*sin(y(2)-y(1));
  dy_dt = f+exch;
end