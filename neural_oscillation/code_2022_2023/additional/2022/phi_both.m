y0 = [ 0.0; 0.0];
opts = odeset('RelTol',2e-13,'AbsTol',1e-100, 'MaxStep',0.1);

d = 0;
alpha = 0;
[T,Y] = ode45(@(t,y)eqn(t,y,d, alpha),[0, 600], y0,opts);
Y = mod(Y, 2*pi); 
plot(T, Y(:,1),'-b',T, Y(:,2),'-g');
ylim([0 2*pi]); yticks([0 pi 2*pi]); yticklabels(["0" "\pi" "2\pi"]);

title(sprintf('d = %.3f', d));

legend({'n = 3, \gamma = 1.01', 'n = 3, \gamma = 1.02'},'Location','northwest');
xlabel("t");
ylabel('\phi', 'Interpreter','tex');

function dy_dt = eqn(t,y,d, alpha)
  n = [3; 3];
  g = [ 1.01; 1.02];
  f = g-sin(y./n);
  exch = [d * sin(y(2)-y(1) - alpha); d * (sin(y(1) - y(2) - alpha))];
  dy_dt = f+exch;
%   x = y./n;
%   dy_dt = g-sin(x);
end
