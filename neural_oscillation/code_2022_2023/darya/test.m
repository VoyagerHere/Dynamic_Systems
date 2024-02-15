%задание параметров
d=0.15;
sigma = 1;

n1 = 1;
n2 = 1;

y0 = [pi/2; pi/2]; %временные рамки для построения
a=10000;
b=14000;

opts = odeset('RelTol',2e-14,'AbsTol',1e-150); %решение системы ДУ
F1 = [1 ; 1];
[T1,Y1] = ode45(@(t,y)eqn(t,y,F1,d,n1,n2), [a, b], y0, opts);

F2 = [0 ; 1];
[T2,Y2] = ode45(@(t,y)eqn(t,y,F2,d,n1,n2), [a, b], y0, opts);

F3 = [1 ; 0];
[T3,Y3] = ode45(@(t,y)eqn(t,y,F3,d,n1,n2), [a, b], y0, opts);

F4 = [0 ; 0];
[T4,Y4] = ode45(@(t,y)eqn(t,y,F4,d,n1,n2), [a, b], y0, opts);

%Y1 = mod(Y1, 2*pi);
%Y2 = mod(Y2, 2*pi);
%Y3 = mod(Y3, 2*pi);
%Y4 = mod(Y4, 2*pi);

%блок проверки решений на соответствие принадлежности к D_i

for(i = 1:length(T1))
if((Y1(i,1) < (pi/2) + sigma && Y1(i,1) > (pi/2) - sigma) || (Y1(i,2) < (pi/2) + sigma && Y1(i,2) > (pi/2) - sigma))
Y1(i,1) = -1;
Y1(i,2) = -1;
end
end

for(i = 1:length(T2))
if((Y2(i,1) > (pi/2) + sigma || Y2(i,1) < (pi/2) - sigma) || (Y2(i,2) < (pi/2) + sigma && Y2(i,2) > (pi/2) - sigma))
Y2(i,1) = -1;
Y2(i,2) = -1;
end
end

for(i = 1:length(T3))
if((Y3(i,1) < (pi/2) + sigma && Y3(i,1) > (pi/2) - sigma) || (Y3(i,2) > (pi/2) + sigma || Y3(i,2) < (pi/2) - sigma))
Y3(i,1) = -1;
Y3(i,2) = -1;
end
end

for(i = 1:length(T4))
if((Y4(i,1) > (pi/2) + sigma || Y4(i,1) < (pi/2) - sigma) || (Y4(i,2) > (pi/2) + sigma || Y4(i,2) < (pi/2) - sigma))
Y4(i,1) = -1;
Y4(i,2) = -1;
end
end



%внесение всех подходящих решений в общий массив

Y = [];
T = [];

Nn = length(T1)+length(T2)+length(T3)+length(T4);
counter_1 = 1;
counter_2 = 1;
counter_3 = 1;
counter_4 = 1;
Y_counter = 1;
minus_one_flag = false;
while(true)
if(counter_1 >= length(T1) || counter_2 >= length(T2) || counter_3 >= length(T3) || counter_4 >= length(T4))
break;
end
curr = [T1(counter_1); T2(counter_2); T3(counter_3); T4(counter_4)];
num = find(abs(curr(:)-min(curr)) < 0.00000005,1);
if(num == 1)
counter_1 = counter_1 + 1;
if (Y1(counter_1,1) > -1)
Y(Y_counter,1) = Y1(counter_1,1);
Y(Y_counter,2) = Y1(counter_1,2);
T(Y_counter) = T1(counter_1);
Y_counter = Y_counter + 1;
continue;
end
end
if(num == 2)
counter_2 = counter_2 + 1;
if (Y2(counter_2,1) > -1)
Y(Y_counter,1) = Y2(counter_2,1);
Y(Y_counter,2) = Y2(counter_2,2);
T(Y_counter) = T2(counter_2);
Y_counter = Y_counter + 1;
continue;
end
end
if(num == 3)
counter_3 = counter_3 + 1;
if (Y3(counter_3,1) > -1)
Y(Y_counter,1) = Y3(counter_3,1);
Y(Y_counter,2) = Y3(counter_3,2);
T(Y_counter) = T3(counter_3);
Y_counter = Y_counter + 1;
continue;
end
end
if(num == 4)
counter_4 = counter_4 + 1;
if (Y4(counter_4,1) > -1)
Y(Y_counter,1) = Y4(counter_4,1);
Y(Y_counter,2) = Y4(counter_4,2);
T(Y_counter) = T4(counter_4);
Y_counter = Y_counter + 1;
continue;
end
end
end
%{
i = 1;
skipFlag = false;
while(true)
if(i >= length(T)-2)
break;
end

if(Y(i+1,1) <= 1.06)
if(Y(i,1) > 1.07)
Y(i+1,1) = Y(i,1);
end
if(Y(i+2,1) > 1.07)
Y(i+1,1) = Y(i+2,1);
end
end

if(Y(i+1,2) <= 1.06)
if(Y(i,2) > 1.07)
Y(i+1,2) = Y(i,2);
end
if(Y(i+2,2) > 1.07)
Y(i+1,2) = Y(i+2,2);
end
end

i = i + 1;
end
%}

newY1 = findpeaks(Y(:,1));
newY2 = findpeaks(Y(:,2));
newT1 = [];
for(i = 1:length(newY1))
newT1(i) = T(find(abs(newY1(i) - Y(:,1)) < 0.000005,1));
end

newT2 = [];
for(i = 1:length(newY2))
newT2(i) = T(find(abs(newY2(i) - Y(:,2)) < 0.000005,1));
end

newY1 = mod(newY1, 2*pi);
newY2 = mod(newY2, 2*pi);

%Y = mod(Y, 2*pi);
%coun = 1;

%while(coun < length(T)-2)
%if(Y(coun,1) > Y(coun+1,1) && abs(Y(coun,1) - Y(coun+1,1)) > 0.05)
%if(Y(coun+2,1) > pi/2)
%Y(coun+1,1) = Y(coun+2,1);
%else
%Y(coun,1) = Y(coun+1,1);
%end
%end
%if(Y(coun,2) > Y(coun+1,2) && abs(Y(coun,2) - Y(coun+1,2)) > 0.05)
%if(Y(coun+2,2) > pi/2)
%Y(coun+1,2) = Y(coun+2,2);
%else
%Y(coun,2) = Y(coun+1,2);
%end
%end
%coun = coun + 1;
%end%}
%построение графика
%disp(length(T))

plot(newT1, newY1,'-b',newT2, newY2,'-r');
ylim([0 2*pi]); yticks([0 pi 2*pi]); yticklabels(["0" "\pi" "2\pi"]);
title('d = 0.15', 'sigma = 1');
legend('n_1=1','n_2=1', 'Interpreter','tex');
xlabel("t");
ylabel('\theta', 'Interpreter','tex');


%задание ДУ

function dy_dt = eqn(~,y,F,d,n1,n2)
g = [1.01; 1.02];
no = [n1; n2];
f = g - sin(y./no);
exch = d*F([2;1]);
dy_dt = f-exch;
end