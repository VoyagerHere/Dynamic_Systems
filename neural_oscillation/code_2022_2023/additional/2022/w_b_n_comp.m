format long g

d = 0;
warning('off');
global n;
global g;
W1 = [];



z = 3:20;
% m = [1.01, 1.02, 1.05, 1.1, 1.15, 1.2];
g = 1.01;


for i = 1:length(z)
    n = z(i);
    disp(g);
    disp(n);
    w1 = f1();
    w2 = arrayfun(@(x)f2(x),d);
    W1(i, 1) = w1;
    W1(i, 2) = w2;
    W1(i, 3) = w3;
end
plot(z, W1(:,1), z, W1(:,2), z, W1(:,3));

% hold on;
% plot(z, W1(:,1));
% plot(z,1./(7*z));
% hold off;

title('\gamma = 1.01 d = 0', 'Interpreter','tex');
% legend('Аналитич. график', 'w_b = 1/(7n)', 'Interpreter','tex');
xlabel('n', 'Interpreter','tex');
ylabel('w', 'Interpreter','tex');

function w1 = f1()
  global n;
  global g;
  T_All = T_Calc(n);
  w1 = (2*pi)/(sum(T_All));
  sum(T_All)
end

function w2 = f2(d)
  global n;
  global g;
  y0 = [0; 0];
  a=0;
  b=10000;
  opts = odeset('RelTol',2e-13,'AbsTol',1e-100);
  [T,Y] = ode45(@(t,y)eqn(t,y,n), [a, b], y0, opts);
  
  TT1 = ms(Y,T,n,1);
  TT1 = sum(TT1);
  w2 = 2*pi/TT1;
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
