sigma = 1; %проходит по значениям
d = 0.0:0.01:0.3; %фиксирован в 0.0015
n = 31;

n1 = 1;
n2 = 1;
k=-500;

y0 = [(0); (0)]; %временные рамки для построения
%y0 = [pi; pi];
a=10000;
b=14000;

w1 = [];
w2 = [];
parfor j = 1:n

%I1=1/(1+exp^(k*(cos(sigma)-sin(y.))));
%I2=1/(1+exp^(k*(cos(sigma)-sin(y.))));

opts = odeset('RelTol',2e-13,'AbsTol',1e-100); %решение системы ДУ
[T,Y] = ode45(@(t,y)eqn(t,y,k,d(j),sigma,n1,n2), [a, b], y0, opts);



NForSearchPiLast1 = find(abs(Y(:,1)-21*pi) < 0.04,1);
NForSearchPiLast2 = find(abs(Y(:,2)-21*pi) < 0.04,1);

TForPiLast1 = T(NForSearchPiLast1);
TForPiLast2 = T(NForSearchPiLast2);


w1 = [w1; vpa((2*pi*10)/(TForPiLast1))];
w2 = [w2; vpa((2*pi*10)/(TForPiLast2))];
end
close all
disp(length( w1))
disp(length( w2))
disp(length(sigma))
plot(d,w1,'r',d,w2,'b')

%построение графика временной реализации

%Y=mod(Y, 2*pi) ;
%plot(T, Y(:,1),'-b', T, Y(:,2),'-r');
%ylim([0 2*pi]); yticks([0 pi 2*pi]); yticklabels(["0" "\pi" "2\pi"]);

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
title('d = 0.15', 'sigma = 1');
legend('n_1=1','n_2=1', 'Interpreter','tex');
xlabel("t");
ylabel('\phi', 'Interpreter','tex');

%задание ДУ

function dy_dt = eqn(~,y,k ,sigma ,d,n1,n2)
g = [1.01; 1.02];
no = [n1; n2];
f = g - cos(y./no);
exch = d*(1./(1+exp(k*(cos(sigma)-sin(y([2;1]))))));
dy_dt = f-exch;
end