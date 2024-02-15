% 
% LOW ACCURATE
% 
warning('off')
y0 = [ 0.0; 0.0];

OUT = []; %(per_num, g_num + n_num)
OUT_ALL = []; %(
global g;
gg = 1.005:0.005:1.2;

nn = [3,5,7,9];
n = 4;
    for m = 1:length(gg)
        g = gg(m);
        opts = odeset('RelTol',2e-13,'AbsTol',1e-100, 'MaxStep',0.1);
        [T,Y] = ode45(@(t,y)eqn(t,y,n),[0, 6000], y0,opts);
        tilda = tld(n);        
        for k = 1:n 
            mid1 = tilda+2*pi*(k-1);
            mid2 = mid1 + 2*pi;
            A1 = find(abs(Y(:,1) - mid1) < 0.04);
            A2 = find(abs(Y(:,1) - mid2) < 0.04);
            AA1 = mean(f4(A1,T));
            AA2 = mean(f4(A2,T));
            OUT(k, m) = AA2 - AA1;
        end
        OUT(:,m) = sort( OUT(:,m),'descend');
    end

hold on;
% plot(gg, OUT(1,1:length(gg)), 'b-')
plot(gg, OUT(2,1:length(gg)), 'r--')
plot(gg, OUT(3,1:length(gg)), 'g-.')
plot(gg, OUT(4,1:length(gg)), 'c--')
% plot(gg, OUT(5,1:length(gg)), 'm-.')
% plot(gg, OUT(6,1:length(gg)), 'b--')
% plot(gg, OUT(7,1:length(gg)), 'r-.')
% plot(gg, OUT(8,1:length(gg)), 'c--')
% plot(gg, OUT(9,1:length(gg)), 'm-.')
hold off;
% legend('T_1', 'T_2', 'T_3', 'T_4', 'T_5', 'T_6', 'T_7', 'T_8', 'T_9', 'Interpreter','tex');
% legend('T_1', 'T_2', 'T_3', 'T_4', 'T_5', 'T_6', 'T_7', 'T_8', 'T_9', 'Interpreter','tex');
legend('T_2', 'T_3', 'T_4', 'Interpreter','tex');
xlabel("g");ylabel("T");
title(sprintf('n = %d', n));

function AA1 = f4(AA1, T)
    AA1 = num2cell(AA1);
    AA1 = T([AA1{:}], :);
end

function tilda = tld(n)
    tilda = n*pi/2 - floor(n/4)*2*pi;
    if abs(tilda - pi/2) < 0.02
        tilda = 3*pi/2;
    else
        tilda = abs(tilda - pi);
    end
end


function dy_dt = eqn(t,y,n)
  global g;
  dy_dt = g-sin(y./n);
end
