gamma = 1.01; n=1;
y0 = [ 0.0; 0.0];

a = 0;
b = 600;
opts = odeset('RelTol',2e-13,'AbsTol',1e-100, 'Refine', 100);
[T,Y] = ode45(@(t,y)eqn(t,y,n), [a, b], y0, opts);
YY = mod(Y, 2*pi); 

% % DELETE < 2pi
END = find(abs(YY(:, 1) - 2*pi) < 0.001);
END = Near_Data(END);
RES_Y = Y(1:END(1,1), 1);

% Y = mod(Y, 2*pi); 
% plot(T, Y(:,1),'-b',T, Y(:,2),'-g');
%  
% Calculate interpolated phi
for k = 1:length(END)-1
  [tau_pr, tau] = find_cross(Y,T,n,k);
  if tau == 0
    tau = tau_pr;
    tau_pr = 0;
  end
  phi = 2*pi * k + (T(END(k,1) + 1:END(k+1,1), 1) - tau_pr)/(tau - tau_pr);
  RES_Y(END(k,1)+1:END(k+1,1), 1) = mod(phi, 2*pi);
end


tilda = tld(n);     
hold on;
yline(tilda,'-', ' $$\tilde{\phi}$$', 'Interpreter', 'LaTeX', 'fontsize', ...
    14, 'LabelHorizontalAlignment', 'right', 'LineWidth', 2);

ylim([0 2*pi]); yticks([0 pi 2*pi]); yticklabels(["0" "\pi" "2\pi"]);
plot(T(1:length(RES_Y),1), RES_Y(:,1),'-b');
grid on; grid minor; 
title(sprintf('n = %d', n));
xlabel("t");
ylabel('\phi', 'Interpreter','tex');
hold off;



function [tau_pr, tau] = find_cross(Y, T, n, k)
  [tau_pr, tau] = ms(Y,T,n,1, k);
end


function OUT = ms_freq(Y,T,n,ind)
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
function [AA1, AA2] = ms(Y,T,n,ind, k)
    tilda = tld(n);
    OUT = [];

    
    mid1 = tilda+2*pi*(k-1);
    mid2 = mid1 + 2*pi;
    
    A1 = find(abs(Y(:,ind) - mid1) < 0.005);
    A2 = find(abs(Y(:,ind) - mid2) < 0.005);
    AA1 = mean(id(A1,T));
    AA2 = mean(id(A2,T));    

%     if phi = [0, 2pi]
    AA2(isnan(AA2))=0; 
end

function AA1 = Near_Data(A1)
  AA1 = [];
  fin = A1(1);
  for k = 1:length(A1)
    if (k ~= length(A1))
      if ((A1(k+1) - A1(k)) < 10)
        fin = A1(k+1);
      else
        AA1(end+1, 1) = fin;
        fin = A1(k+1);
      end
    else
      AA1(end+1,1) = A1(k);
    end
  end
end


function tilda = tld(n)
    tilda = n*pi/2 - floor(n/4)*2*pi;
    if tilda == pi/2
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
  g = 1.01;
  dy_dt = g-sin(y./n);
end