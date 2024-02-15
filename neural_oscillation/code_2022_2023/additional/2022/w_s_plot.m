
warning('off')
d = 0;
% d = 3;
% d = 0:0.2:3;

global epsilon; 
epsilon = 0.002;
global n;
global g;
W1 = [];
W2 = [];


z = 3:10;
m = [1.01, 1.02, 1.05, 1.1, 1.15, 1.2];
for k = 1:length(m)
    g = m(k);
    for i = 1:length(z)
        n = z(i);
        disp(g);
        disp(n);
        w1 = arrayfun(@(x)f1(x),d);
        W1(i, k) = w1;
    end
end

plot(z, W1(:,1), z, W1(:,2), z, W1(:,3), z, W1(:,4), z, W1(:,5), z, W1(:,6));
legend('\gamma = 1.01', '\gamma = 1.02','\gamma = 1.05', '\gamma = 1.1', '\gamma = 1.15','\gamma = 1.2','Interpreter','tex');
xlabel('n', 'Interpreter','tex');
ylabel('w_s', 'Interpreter','tex');



function w1 = f1(d)
  global  n;
  i = 0;
  y0 = [0; 0];
  a=0;
  b=10000;
  opts = odeset('RelTol',2e-13,'AbsTol',1e-100, 'Refine', 100);
  [T,Y] = ode45(@(t,y)eqn(t,y,d), [a, b], y0, opts);
  TT1 = ms(Y,T,n,1);
  if (size(TT1,1) < 4)
      i = 1;
  end
  TT1 = sort(TT1,'descend');
  TT1 = TT1(2:size(TT1));
  TT1 = mean(TT1);
  w1 = 2*pi/TT1;
end

function dy_dt = eqn(t,y,n, g)
  global g;
  global n;
  dy_dt = g-sin(y./n);
end

function OUT = ms(Y,T,n,ind)
    tilda = tld(n);
    OUT = [];
    error = 0.01;

    for k = 1:n 
        mid1 = tilda+2*pi*(k-1);
        mid2 = mid1 + 2*pi;
        
        A1 = find(abs(Y(:,ind) - mid1) < error);
        A2 = find(abs(Y(:,ind) - mid2) < error);
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