
d = 0;
warning('off');
% d = 0.006;

global n;
global g;
W1 = [];



z = 3:10;
m = [1.01, 1.02, 1.05, 1.1, 1.15, 1.2];
for k = 1:length(m)
    g = m(k);
    for i = 1:length(z)
        n = z(i);
        disp(g);
        disp(n);
        w1 = f1();
        W1(i, k) = w1;
    end
end
plot(z, W1(:,1), z, W1(:,2), z, W1(:,3), z, W1(:,4), z, W1(:,5), z, W1(:,6));
legend('\gamma = 1.01', '\gamma = 1.02','\gamma = 1.05', '\gamma = 1.1', '\gamma = 1.15','\gamma = 1.2','Interpreter','tex');
xlabel('n', 'Interpreter','tex');
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

function AA1 = id(AA1, T)
    AA1 = num2cell(AA1);
    AA1 = T([AA1{:}], :);
end

