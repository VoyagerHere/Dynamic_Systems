d = 0:0.001:0.5;
% d = 0.028;
% d = 0:0.2:3;
global epsilon;
epsilon = 0.002;

global begin_error;
begin_error = 0.1;


global D;
D = [];
global n1;
global n2;
n1 = 3;
n2 = 3;

[w1, w2] = arrayfun(@(x)f1(x),d);
d_cr = D(1,1);

hold on
plot(d, w1,'b');
plot(d, w2, 'g');

xline(d_cr,'-', ' d_{cr}^{s}', 'Interpreter','tex', 'LineWidth', 2);
hold off

legend('w_{s}^{1}','w_{s}^{2}', 'Interpreter','tex');
title('n_1 = n_2 = 3', 'Interpreter','tex');
xlabel("d");
ylabel('w_s', 'Interpreter','tex');

ylim([1 1.4]);

function [w1,w2] = f1(d)
  global begin_error;
  global  n1;
  global  n2;
  global D;
  global epsilon;
  y0 = [0; 0];
  a=0;
  b=600;
  opts = odeset('RelTol',2e-13,'AbsTol',1e-100, 'Refine', 100);
  [T,Y] = ode45(@(t,y)eqn(t,y,d), [a, b], y0, opts);

  TT2 = ms(Y,T,n2,2);
  TT1 = ms(Y,T,n1,1);
  TT2 = sort(TT2,'descend');
  TT1 = sort(TT1,'descend');

  TT2 = TT2(2:size(TT2));
  TT1 = TT1(2:size(TT1));
  TT2 = mean(TT2);
  TT1 = mean(TT1);
  
  

  w2 = vpa(2*pi/TT2);
  w1 = vpa(2*pi/TT1);
  if (abs(w2 - w1) < epsilon)
    if (d > begin_error)
      D(end+1,:) = d;
    end
  end
end


function dy_dt = eqn(~,y,d)
  global  n1;
  global  n2;
  g = [1.01; 1.02];
  no = [n1; n2];
  f = g-sin(y./no);
  exch = [d;-d]*sin(y(2)-y(1));
  dy_dt = f+exch;
end


function OUT = ms(Y,T,n,ind)
    tilda = tld(n);
    OUT = [];

    for k = 1:n
        mid1 = tilda+2*pi*(k-1);
        mid2 = mid1 + 2*pi;
        
        A1 = find(abs(Y(:,ind) - mid1) < 0.002);
        A2 = find(abs(Y(:,ind) - mid2) < 0.002);
        AA1 = mean(id(A1,T));
        AA2 = mean(id(A2,T));
        OUT(end+1, :) = AA2 - AA1;
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

function AA1 = id(AA1, T)
    AA1 = num2cell(AA1);
    AA1 = T([AA1{:}], :);
end