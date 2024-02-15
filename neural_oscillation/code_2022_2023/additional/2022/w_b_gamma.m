d = 0;
warning('off');

global n;
global g;
W1 = [];


gg = linspace(1.01, 2.01, 10);
z = [3; 6; 9];
for k = 1:length(z)
    n = z(k);
    for i = 1:length(gg)
        g = gg(i);
        disp(g);
        disp(n);
        w1 = f1();
        W1(i, k) = w1;
    end
end
plot(gg, W1(:,1),gg, W1(:,2), gg, W1(:,3));
legend('n=3','n=6', 'n=9');
xlabel('\gamma', 'Interpreter','tex');
ylabel('w_b', 'Interpreter','tex');

function w1 = f1()
  global n;
  global g;
  T_All = T_Calc(n);
  w1 = (2*pi)/(sum(T_All));
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

function tilda = tld(n)
    tilda = n*pi/2 - floor(n/4)*2*pi;
    if abs(tilda - pi/2) < 0.02
        tilda = 3*pi/2;
    else
        tilda = abs(tilda - pi);
    end
end
