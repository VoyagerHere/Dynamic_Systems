d = 0:0.001:0.25;
% d = 6;
warning('off')


global epsilon;
% epsilon = 0.002;
epsilon = 0.0001;



global D;
D = [];


global n1;
global n2;
nn1 = 2:10;

D_CR = zeros(length(nn1), 1);

for k = 1:length(nn1)
    n1 = nn1(k);
    n2 = n1;
    
    D = [];
    disp(n1);

    [w1, w2] = arrayfun(@(x)f1(x),d);
    d_cr = find_cr(D);
    D_CR(k, 1) = d_cr;
end

plot(nn1, D_CR);
xlabel('(n_1, n_2)', 'Interpreter','tex');
ylabel('d_{cr}^{b}', 'Interpreter','tex');
xticklabels({'(1,1)', '(2,2)', '(3,3)','(4,4)','(5,5)','(6,6)','(7,7)','(8,8)','(9,9)','(10,10)'})

function [w1, w2] = f1(d)
  global D;
  global epsilon;
  global n1;
  global n2;
  y0 = [0; 0];
  a=0;
  b=10000;
  opts = odeset('RelTol',2e-13,'AbsTol',1e-100, 'Refine', 100);
  [T,Y] = ode45(@(t,y)eqn(t,y,d), [a, b], y0, opts);
  TT2 = ms(Y,T,n2,2);
  TT1 = ms(Y,T,n1,1);

  TT2 = sum(TT2);
  TT1 = sum(TT1);
 
  w2 = 2*pi/TT2;
  w1 = 2*pi/TT1;
  disp(d)
  if (abs(w2 - w1) < epsilon)
    len = size(D,1);
    D(len+1,1) = d;
    D(len+1,2) = w1;
    return;
  end
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





