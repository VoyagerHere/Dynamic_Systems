y0 = [ 0.0; 0.0; 0.0];


g1 = 1.01;
% d value
d1 = 0.08;
% d2 = 0.02;
d2 = d1;
delta = 0.01;

% alpha value
alpha = 0;


g = [g1, g1 + delta, g1 + 2*delta];
opts = odeset('RelTol',2e-13,'AbsTol',1e-100, 'MaxStep',0.1);
[T,Y] = ode45(@(t,y)eqn(t,y,g, d1,d2,alpha),[2000, 2500], y0,opts);
Y = mod(Y, 2*pi);
plot(T, Y(:,1),'-b',T, Y(:,2),'-g', T, Y(:,3),'-m');
ylim([0 2*pi]); yticks([0 pi 2*pi]); yticklabels(["0" "\pi" "2\pi"]);

title(sprintf('Δ = %g, d_1 = %g, d_2 = %g', delta, d1, d2));

legend({'\phi_1, n_1 = 3, \gamma_1 = 1.01', '\phi_2, n_2 = 3, \gamma_2 = 1.02', '\phi_3, n_3 = 3, \gamma_3 = 1.03'},'Location','northeast','FontSize', 10,'FontWeight','bold');
xlabel("t",'FontSize', 12,'FontWeight','bold');
ylabel('\phi', 'Interpreter','tex','FontSize', 15,'FontWeight','bold');


function dy_dt = eqn(t,y,g, d1,d2, alpha)
  n = [3; 3; 3];
  g = [g(1); g(2); g(3)];
  f = g-sin(y./n);
  exch = [d1*sin(y(2) - y(1) - alpha); d1*sin(y(1) - y(2) - alpha) + d2*sin(y(3) - y(2) - alpha); d2*sin(y(2) - y(3) - alpha)];
  dy_dt = f+exch;
end
