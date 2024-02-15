format long g

d = 0;
warning('off');
global n;
global g;
W1 = [];



z = 3:10;
% m = [1.01, 1.02, 1.05, 1.1, 1.15, 1.2];
g = 1.01;


for i = 1:length(z)
    n = z(i);
    disp(g);
    disp(n);
    w1 = f1();
    w2 = f2();
    w3 = f3();
    W1(i, 1) = w1;
    W1(i, 2) = w2;
    W1(i, 3) = w3;
end
plot(z, W1(:,1), z, W1(:,2), z, W1(:,3));
title('\gamma = 1.01 d = 0', 'Interpreter','tex');
legend('w_b', 'w_c', 'w_i', 'Interpreter','tex');
xlabel('n', 'Interpreter','tex');
ylabel('w', 'Interpreter','tex');
xlim([3 10])

function w1 = f1()
% burst
  global n;
  global g;
  T_All = T_Calc(n);
  w1 = (2*pi)/(sum(T_All));
  sum(T_All)
end


function w1 = f2()
  global  n;
  y0 = [0; 0];
  a=0;
  b=10000;
  opts = odeset('RelTol',2e-13,'AbsTol',1e-100, 'Refine', 100);
  [T,Y] = ode45(@(t,y)eqn(t,y,n), [a, b], y0, opts);
  w1 = vpa((Y(end,1) - Y(502,1))/(b*2*pi));
end


function w1 = f3()
  global  n;
  y0 = [0; 0];
  a=0;
  b=100000;
  opts = odeset('RelTol',2e-13,'AbsTol',1e-100, 'Refine', 100);
  [T,Y] = ode45(@(t,y)eqn(t,y,n), [a, b], y0, opts);
  num_1 = ms_2(Y,n,1);
  w1 = vpa(2*pi*(num_1)/b);
end

function num = ms_2(Y, n,ind)
    tilda = tld(n);
    OUT = [];
    k = 1;
    while 1
        mid1 = tilda+2*pi*(k-1);
        A1 = find(abs(Y(:,ind) - mid1) < 0.02);
        if isempty(A1)
          break;
        end
        k = k + 1;
    end
    num = k - 1;
end

function OUT = T_Calc(n)
  global g;
  OUT = [];
  tilda = tld(n);
  for k = 1:n 
      phi1 = tilda+2*pi*(k-1);
      phi2 = phi1 + 2*pi;
      
      x1 =  (2 * n * (atan((g*(tan(phi1/(2*n))- 1)/sqrt(g^2 - 1)))))/sqrt(g^2 - 1);
      x2 =  (2 * n * (atan((g*(tan(phi2/(2*n))- 1)/sqrt(g^2 - 1)))))/sqrt(g^2 - 1);
      T = x2 - x1;
      x3 =  (2 * n * (atan((g*(tan(phi2/(2*n))- 1)/sqrt(g^2 - 1)))+pi))/sqrt(g^2 - 1);
      TT = x3 - x1;
      if (T > 0)
        OUT(end+1, 1) = T;
      else
        OUT(end+1, 1) = TT;
      end
  end
end





function OUT = ms(Y,T,n,ind)
    tilda = tld(n);
    OUT = [];
%     error = 0.002;
    error = 0.05;
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




function dy_dt = eqn(t,y,n)
  global g;
  dy_dt = g-sin(y./n);
end



