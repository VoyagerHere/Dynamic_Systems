sigma=1;
d = 0.0:0.01:15.0;
n = 1501;

n1 = 1;
n2 = 1;


y0 = [(pi/2); (pi/2)]; 
a=10000;
b=12000;

w1 = [];
w2 = [];

parfor j = 1:n


opts = odeset('RelTol',2e-13,'AbsTol',1e-100); %решение системы ДУ
F = [1 ; 1];
[T1,Y1] = ode45(@(t,y)eqn(t,y,F,d,n1,n2), [a, b], y0, opts);

F = [0 ; 1];
[T2,Y2] = ode45(@(t,y)eqn(t,y,F,d,n1,n2), [a, b], y0, opts);

F = [1 ; 0];
[T3,Y3] = ode45(@(t,y)eqn(t,y,F,d,n1,n2), [a, b], y0, opts);

F = [0 ; 0];
[T4,Y4] = ode45(@(t,y)eqn(t,y,F,d,n1,n2), [a, b], y0, opts);



N = min([length(T1),length(T2),length(T3),length(T4)]);
for(i = 1:N)
if((Y1(i,1) < (pi/2) + sigma(j) && Y1(i,1) > (pi/2) - sigma(j)) || (Y1(i,2) < (pi/2) + sigma(j) && Y1(i,2) > (pi/2) -sigma(j)))
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





Y = [];

for(i = 1:N)
if(Y1(i,1) > -1)
Y(i,1) = Y1(i,1);
Y(i,2) = Y1(i,2);
end
if(Y2(i,1) > -1)
Y(i,1) = Y2(i,1);
Y(i,2) = Y2(i,2);
end
if(Y3(i,1) > -1)
Y(i,1) = Y3(i,1);
Y(i,2) = Y3(i,2);
end
if(Y4(i,1) > -1)
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



NForSearchPiLast1 = find(abs(Y(:,1)-21*pi) < 0.004,1);
NForSearchPiLast2 = find(abs(Y(:,2)-21*pi) < 0.004,1);

TForPiLast1 = T(NForSearchPiLast1);
TForPiLast2 = T(NForSearchPiLast2);


w1 = [w1; vpa((2*pi*10)/(TForPiLast1))];
w2 = [w2; vpa((2*pi*10)/(TForPiLast2))];
end

close all
disp(length( w1))
disp(length( w2))

plot(d(1:length( w1)),w1,'r',d(1:length( w2)),w2,'b')


title('d = 0.15', 'sigma = 1');
legend('n_1=1','n_2=1', 'Interpreter','tex');
xlabel("t");
ylabel('\phi', 'Interpreter','tex');



function dy_dt = eqn(~,y ,d,n1,n2)
g = [1.01; 1.02];
no = [n1; n2];
f = g - sin(y./no);
exch = d*F([2;1]);
dy_dt = f-exch;
end