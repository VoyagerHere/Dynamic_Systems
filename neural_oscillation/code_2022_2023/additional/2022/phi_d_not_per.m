y0 = [ 0.0; 0.0 ];
opts = odeset('RelTol',2e-13,'AbsTol',1e-100,'Refine', 100);
a = 0;
b = 10000;


d = 0;
[T,Y] = ode45(@(t,y)eqn(t,y,d), [a, b], y0, opts);

delta1 = Y(end,1) - Y(500,1);
delta2 = Y(end,2) - Y(500,2);
w1 = vpa((Y(end,1) - Y(500,1))/(b*2*pi));
w2 = vpa((Y(end,2) - Y(500,2))/(b*2*pi));


% % PLOT

hold on;

% 
% title('\gamma = 1.01, d = 0','Interpreter','tex');
% Y = mod(Y, 2*pi); 
% ylim([0 2*pi]); yticks([0 pi 2*pi]); yticklabels(["0" "\pi" "2\pi"]);

plot(T, Y(:,1),'-g',T, Y(:,2),'-b');

% 
% ellipse_params = fit_ellipse([x1, x1], [0, y1], gca);
% ellipse(ellipse_params.a, ellipse_params.b, ellipse_params.phi, ...
%     ellipse_params.X0, ellipse_params.Y0, 'r');


hold off;


title('d = 0, n_1 = 3, n_2 = 3','Interpreter','tex');
legend({'$n_1 = 3, \gamma_1 = 1.01, \Omega_1 = 0.1417$', '$n_2 = 3, \gamma_2 = 1.02, \Omega_2 = 0.2009$'}, 'Interpreter', 'latex', 'Location', 'northwest');
xlabel("t");
ylabel('\phi', 'Interpreter','tex');

% % \ PLOT

function dy_dt = eqn(~,y,d)
  n1 = 3;
  n2 = 3;
  n = [n1; n2];
  g = [ 1.01; 1.02];
  f = g-sin(y./n);
  exch = [d;-d]*sin(y(2)-y(1));
  dy_dt = f+exch;
end
