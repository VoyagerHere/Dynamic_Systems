gamma = 1.01; n=1;
y0 = [ 0.0; 0.0];

opts = odeset('RelTol',2e-13,'AbsTol',1e-100, 'MaxStep',0.1);


[T,Y] = ode45(@(t,x)gamma-sin(x/n), [0,600], y0, opts);
YY = mod(Y, 2*pi); 


% % DELETE < 2pi
END = find(abs(YY(:, 1) - 2*pi) < 0.005);
RES_Y = Y(1:END(1,1), 1);

 
% Calculate interpolated phi
for k = 1:length(END)-1
  [tau_pr, tau] = find_cross(Y,T,n,k);
  if tau == 0
    tau = tau_pr;
    tau_pr = 0;
  end
  phi = 2*pi * k + (T(END(k,1) + 1:END(k+1,1), 1) - tau_pr)/(tau - tau_pr);
  RES_Y(END(k,1)+1:END(k+1,1), 1) = phi;
end


tilda = tld(n);     
RES_Y = mod(RES_Y, 2*pi); 
hold on;
yline(tilda,'-', ' $$\tilde{\phi}$$', 'Interpreter', 'LaTeX', 'fontsize', ...
    14, 'LabelHorizontalAlignment', 'right', 'LineWidth', 2);

ylim([0 2*pi]); yticks([0 pi 2*pi]); yticklabels(["0" "\pi" "2\pi"]);
plot(T, Y(:,1),'-b');
grid on; grid minor; 
title(sprintf('n = %d', n));
xlabel("t");
ylabel('\phi', 'Interpreter','tex');
hold off;



function [tau_pr, tau] = find_cross(Y, T, n, k)
  [tau_pr, tau] = ms(Y,T,n,1, k);
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