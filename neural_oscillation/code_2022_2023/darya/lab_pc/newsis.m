d=0.15;
sigma = 1;

n1 = 1;
n2 = 1;
k=-500;

y0 = [(pi/2); (pi/2)]; %временные рамки для построения
%y0 = [pi; pi];
a=2000;
b=3000;

%I1=1/(1+exp^(k*(cos(sigma)-sin(y.))));
%I2=1/(1+exp^(k*(cos(sigma)-sin(y.))));

opts = odeset('RelTol',2e-13,'AbsTol',1e-100); %решение системы ДУ
[T,Y] = ode45(@(t,y)eqn(t,y,k,sigma,d,n1,n2), [a, b], y0, opts);


%построение графика временной реализации
close all
%Y=mod(Y, 2*pi) ;
plot(T, Y(:,1),'-b', T, Y(:,2),'-r');
ylim([0 2*pi]); yticks([0 pi 2*pi]); yticklabels(["0" "\pi" "2\pi"]);

title('d = 0.15', 'sigma = 1');
legend('n_1=1','n_2=1', 'Interpreter','tex');
xlabel("t");
ylabel('\theta', 'Interpreter','tex');
%ylabel('\phi', 'Interpreter','tex');

%построение ФП
%plot(Y(:,1), Y(:,2),'-r');
%ylim([0 2*pi]); yticks([0 pi 2*pi]); yticklabels(["0" "\pi" "2\pi"]);
%xlim([0 2*pi]); xticks([0 pi 2*pi]); xticklabels(["0" "\pi" "2\pi"]);
%title('d=0.15', 'sigma = 1');
%legend('n_1=1','n_2=1', 'Interpreter','tex');
%xlabel('\theta_1', 'Interpreter','tex');
%ylabel('\theta_2', 'Interpreter','tex');
%xlabel('\phi_1', 'Interpreter','tex');
%ylabel('\phi_2', 'Interpreter','tex');

%задание ДУ

function dy_dt = eqn(~,y,k ,sigma ,d,n1,n2)
g = [1.01; 1.02];
no = [n1; n2];
f = g - sin(y./no);
exch = d*(1./(1+exp(k*(cos(sigma)-sin(y([2;1]))))));
dy_dt = f-exch;
end