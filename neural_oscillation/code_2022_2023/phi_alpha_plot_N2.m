y0 = [ 0.0; 0.0];



% d value
d = 0.08;
delta = 0.01;

% alpha value
alpha = 2*pi/3;
alpha_text = '2π/3';
n = [3; 3];

el1_text = '\phi_1, n_1 = 3, \gamma_1 = 1.01';
el2_text = '\phi_2, n_2 = 3, \gamma_2 = 1.02';

g = [g1; g1 + delta];
opts = odeset('RelTol',2e-13,'AbsTol',1e-14, 'MaxStep',0.1);
[T,Y] = ode45(@(t,y)eqn(g, n, t, y, d,alpha),[2000, 4000], y0,opts);
Y = mod(Y, 2*pi);
plot(T, Y(:,1),'-b',T, Y(:,2),'-g');
ylim([0 2*pi]); yticks([0 pi 2*pi]); yticklabels(["0" "\pi" "2\pi"]);

title(sprintf('Δ = %g, d = %g, α = %s', delta, d, alpha_text));

legend({el1_text, el2_text},'Location','northeast','FontSize', 10,'FontWeight','bold');
xlabel("t",'FontSize', 12,'FontWeight','bold');
ylabel('\phi', 'Interpreter','tex','FontSize', 15,'FontWeight','bold');


function dy_dt = eqn(g, n, ~, y, d, alpha)
    f = g - sin(y ./ n);
    exch = [d * sin(y(2) - y(1) - alpha); d * (sin(y(1) - y(2) - alpha))];
    dy_dt = f + exch;
end
