d = 0:0.001:0.9;
% d = 3;
% d = 0:0.2:3;

global d_cr;
global epsilon; 
epsilon = 0.002;
global begin_error;
begin_error = 0.1;
global D;
D = [];
global n;
global DCR;
DCR = [];



z = 2:10;
for n=z 
    w1 = arrayfun(@(x)f1(x),d);
%     w2 = arrayfun(@(l)f2(l),d);
    d_cr = D(1,1);
    DCR(end+1,1) = d_cr;
    D = [];
end

plot(z, DCR);
xlabel('(n_1, n_2)', 'Interpreter','tex');
ylabel('d_{cr}^{s}', 'Interpreter','tex');
xticklabels({'(2,2)', '(3,3)','(4,4)','(5,5)','(6,6)','(7,7)','(8,8)','(9,9)','(10,10)'})



function w1 = f1(d)
  global  n;
  global D;
  global epsilon;
  global begin_error;
  i = 0;
  y0 = [0; 0];
  a=0;
  b=10000;
  opts = odeset('RelTol',2e-13,'AbsTol',1e-100, 'Refine', 10);
  [T,Y] = ode45(@(t,y)eqn(t,y,d), [a, b], y0, opts);
  TT2 = ms(Y,T,n,2);
  TT1 = ms(Y,T,n,1);

  if (size(TT1,1) < 4)
      i = 1;
  end
  TT2 = sort(TT2,'descend');
  TT1 = sort(TT1,'descend');

  TT2 = TT2(2:size(TT2));
  TT1 = TT1(2:size(TT1));
  TT2 = mean(TT2);
  TT1 = mean(TT1);
  

  w2 = 2*pi/TT2;
  w1 = 2*pi/TT1;
  if (abs(w2 - w1) < epsilon)
    if (d > begin_error)
      D(end+1,:) = d;
    end
  end
end


function dy_dt = eqn(~,y,d)
  global  n;
  g = [1.01; 1.02];
  no = [n; n];
  f = g-sin(y./no);
  exch = [d;-d]*sin(y(2)-y(1));
  dy_dt = f+exch;
end


function OUT = ms(Y,T,n,ind)
    tilda = tld(n);
    OUT = [];
    error = 0.002;

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