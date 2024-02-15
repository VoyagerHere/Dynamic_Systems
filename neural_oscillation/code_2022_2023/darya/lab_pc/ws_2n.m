sigma=1;
d = 0.0:0.01:15.0;
n = 1501;


%global n3;

%n3 = 1;

w1 = [];
w2 = [];


parfor j = 1:n
y0 = [(pi/2); (pi/2)];
a=10000;
b=12000;
n1 = 1;
n2 = 1;

%PHIFirst = 0;

%opts = odeset('RelTol',2e-13,'AbsTol',1e-100);
%[T,Y] = ode45(@(t,y)eqn(t,y,sigma, d), [a, b], y0, opts);
opts = odeset('RelTol',2e-13,'AbsTol',1e-100); %решение системы ДУ
F = [1 ; 1];
[T1,Y1] = ode45(@(t,y)eqn(t,y,F,d,n1,n2), [a, b], y0, opts);

F = [0 ; 1];
[T2,Y2] = ode45(@(t,y)eqn(t,y,F,d,n1,n2), [a, b], y0, opts);

F = [1 ; 0];
[T3,Y3] = ode45(@(t,y)eqn(t,y,F,d,n1,n2), [a, b], y0, opts);

F = [0 ; 0];
[T4,Y4] = ode45(@(t,y)eqn(t,y,F,d,n1,n2), [a, b], y0, opts);


%блок проверки решений на соответствие принадлежности к D_i
N = min([length(T1),length(T2),length(T3),length(T4)]);
for(i = 1:N)
if((Y1(i,1) < (pi/2) + sigma(j) && Y1(i,1) > (pi/2) - sigma(j)) || (Y1(i,2) < (pi/2) + sigma(j) && Y1(i,2) > (pi/2) - sigma(j)))
Y1(i,1) = -1;
Y1(i,2) = -1;
end

if((Y2(i,1) > (pi/2) + sigma(j) || Y2(i,1) < (pi/2) - sigma(j)) || (Y2(i,2) < (pi/2) + sigma(j) && Y2(i,2) > (pi/2) - sigma(j)))
Y2(i,1) = -1;
Y2(i,2) = -1;
end

if((Y3(i,1) < (pi/2) + sigma(j) && Y3(i,1) > (pi/2) - sigma(j)) || (Y3(i,2) > (pi/2) + sigma(j) || Y3(i,2) < (pi/2) - sigma(j)))
Y3(i,1) = -1;
Y3(i,2) = -1;
end

if((Y4(i,1) > (pi/2) + sigma(j) || Y4(i,1) < (pi/2) - sigma(j)) || (Y4(i,2) > (pi/2) + sigma(j) || Y4(i,2) < (pi/2) - sigma(j)))
Y4(i,1) = -1;
Y4(i,2) = -1;
end
end



%внесение всех подходящих решений в общий массив

Y = [];

for(i = 1:N)
if(Y1(i,1) > 0)
Y(i,1) = Y1(i,1);
Y(i,2) = Y1(i,2);
end
if(Y2(i,1) > 0)
Y(i,1) = Y2(i,1);
Y(i,2) = Y2(i,2);
end
if(Y3(i,1) > 0)
Y(i,1) = Y3(i,1);
Y(i,2) = Y3(i,2);
end
if(Y4(i,1) > 0)
Y(i,1) = Y4(i,1);
Y(i,2) = Y4(i,2);
end
end

T=[];
if(N == length(T1))
T = T1;
end
if(N == length(T2))
T = T2;
end
if(N == length(T3))
T = T3;
end
if(N == length(T4))
T = T4;
end

PHIFirst = 0;


NForSearchPiLast1 = find(abs(Y(:,1)-40*pi) < 0.04,1);
NForSearchPiLast2 = find(abs(Y(:,2)-40*pi) < 0.04,1);

TForPiLast1 = T(NForSearchPiLast1);
TForPiLast2 = T(NForSearchPiLast2);

w1 =[w1, vpa((2*pi*20)/(TForPiLast1))];
w2 =[w2, vpa((2*pi*20)/(TForPiLast2))];

end
%d = 0.0015*ones(1,11);
disp(length(d))
disp(length(w1))
disp(length(w2))
%hold on
%axis equal;
close all;
plot(d, w1,'g',d, w2, 'b'); %строит по д2, если менять, то нужно поменять тут на д1

%plot(sigma, w1,'g');

%w = [];
%for(i = 1:n)
%w = [w;w1(i)/w2(i)];
%end
%plot(sigma,w,'r')
%plot(d, w2, 'b')
%d_cr = D(1,1);

ylim([0 2*pi]); %ограничение по y
%xline(d_cr,'-', ' d_{cr}^{b}', 'Interpreter','tex', 'LineWidth', 2);
%hold off

legend('w_{s}^{1}','w_{s}^{2}', 'Interpreter','tex');
title('n_1 = n_2 = 1', 'Interpreter','tex');
xlabel("d");
ylabel('w_s', 'Interpreter','tex');

function dy_dt = eqn(~,y,F,d,n1,n2)
g = [1.01; 1.02];
no = [n1; n2];
f = g - sin(y./no);
exch = d*F([2;1]);
dy_dt = f-exch;
end