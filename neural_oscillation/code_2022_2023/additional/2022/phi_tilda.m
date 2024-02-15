gamma = 1.01; n=2;

opts = odeset( 'AbsTol',1e-6 ,'RelTol',1e-9,'MaxStep', 0.1 );
opts = odeset(opts,'Events',@(t,y)event(t,y));
sol = ode45(@(t,x)gamma-sin(x/n), [0,600], pi, opts);
tfs = [sol.xe sol.x(end)];
N = length(tfs);
clf;
t0 = 0;

tilda = tld(n);     


for i=1:N
    tf = tfs(i);
    t = linspace(t0+1e-2,tf-1e-2,150);
    y = deval(sol,t);
    y = mod(y,2*pi); 
    plot(t, y, 'k');
    hold on; 
    t0=tf;
end

yline(tilda,'-', ' $$\tilde{\phi}$$', 'Interpreter', 'LaTeX', 'fontsize', ...
    14, 'LabelHorizontalAlignment', 'right', 'LineWidth', 2);

hold off;

ylim([0 2*pi]); yticks([0 pi/2 pi 3*pi/2 2*pi]); yticklabels(["0" "\pi/2" "\pi" "3\pi/2" "2\pi"]); 
grid on; grid minor; 
title(sprintf('n = %d', n));
xlabel("t");
ylabel('\phi', 'Interpreter','tex');


function tilda = tld(n)
    tilda = n*pi/2 - floor(n/4)*2*pi;
    if tilda == pi/2
        tilda = 3*pi/2;
    else
        tilda = abs(tilda - pi);
    end
end

function [ res, term, dir ] = event(t,y)
    y = mod(y+pi,2*pi)-pi;
    res =  y ; 
    dir = 1;
    term = 0;
end
