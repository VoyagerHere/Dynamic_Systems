sigma = 1; %проходит по значениям
d = 0.0:0.01:0.3; %фиксирован в 0.0015
n = 31;

w1 = zeros(n,2);
w2 = zeros(n,2);
fprintf('Progress:\n');
fprintf(['\n' repmat('.',1,n) '\n\n']);

parfor j = 1:n
y0 = [pi/2; pi/2];
a=10000;
b=14000;
n1 = 1;
n2 = 1;

PHIFirst = 0;

%opts = odeset('RelTol',2e-13,'AbsTol',1e-100);
%[T,Y] = ode45(@(t,y)eqn(t,y,sigma, d), [a, b], y0, opts);
opts = odeset('RelTol',2e-13,'AbsTol',1e-150); %решение системы ДУ
F1 = [1 ; 1];
[T1,Y1] = ode45(@(t,y)eqn(t,y,F1,d(j),n1,n2), [a, b], y0, opts);

F2 = [0 ; 1];
[T2,Y2] = ode45(@(t,y)eqn(t,y,F2,d(j),n1,n2), [a, b], y0, opts);

F3 = [1 ; 0];
[T3,Y3] = ode45(@(t,y)eqn(t,y,F3,d(j),n1,n2), [a, b], y0, opts);

F4 = [0 ; 0];
[T4,Y4] = ode45(@(t,y)eqn(t,y,F4,d(j),n1,n2), [a, b], y0, opts);


%блок проверки решений на соответствие принадлежности к D_i

Y1 = mod(Y1, 2*pi);
Y2 = mod(Y2, 2*pi);
Y3 = mod(Y3, 2*pi);
Y4 = mod(Y4, 2*pi);



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
counter_1 = 1;
counter_2 = 1;
counter_3 = 1;
counter_4 = 1;
Y_counter = 1;
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

%Y = mod(Y, 2*pi);
coun = 1;

while(coun < length(T)-2)
if(Y(coun,1) > Y(coun+1,1) && abs(Y(coun,1) - Y(coun+1,1)) > 0.05)
if(Y(coun+2,1) > pi/2)
Y(coun+1,1) = Y(coun+2,1);
else
Y(coun,1) = Y(coun+1,1);
end
end
if(Y(coun,2) > Y(coun+1,2) && abs(Y(coun,2) - Y(coun+1,2)) > 0.05)
if(Y(coun+2,2) > pi/2)
Y(coun+1,2) = Y(coun+2,2);
else
Y(coun,2) = Y(coun+1,2);
end
end
coun = coun + 1;
end


x1 = find(Y(:,1) < 0.05);
x2 = find(Y(:,2) < 0.05);
sort(x1);
sort(x2);
w1_tmp = zeros(1,2);
w2_tmp = zeros(1,2);
if(isempty(x1))
w1_tmp(1) = 0;
w1_tmp(2) = d(j);
else
w1_tmp(1) = vpa((length(x1)*2*pi)/T(x1(end)));
w1_tmp(2) = d(j);
end

if(isempty(x2))
w2_tmp(1) = 0;
w2_tmp(2) = d(j);
else
w2_tmp(1) = vpa((length(x2)*2*pi)/T(x2(end)));
w2_tmp(2) = d(j);
end
w1(j,:) = w1_tmp;
w2(j,:) = w2_tmp;

fprintf('\b|\n');
end
%d = 0.0015*ones(1,19);
%hold on
%axis equal;
close all;
plot(w1(:,2), w1(:,1),'g',w2(:,2), w2(:,1), 'b'); 
%w = [];
%for(i = 1:n)
%w = [w;w1(i)/w2(i)];
%end
%plot(sigma,w,'r')
%plot(d, w2, 'b')
%d_cr = D(1,1);

%ylim([0 2*pi]); %ограничение по y
%xline(d_cr,'-', ' d_{cr}^{b}', 'Interpreter','tex', 'LineWidth', 2);
%hold off

legend('w_{s}^{1}','w_{s}^{2}', 'Interpreter','tex');
title('n_1 = n_2 = 1; sigma = 1', 'Interpreter','tex');
xlabel("d");
ylabel('w_s', 'Interpreter','tex');

function dy_dt = eqn(~,y,F,d,n1,n2)
g = [1.01; 1.02];
no = [n1; n2];
f = g - sin(y./no);
exch = d*F([2;1]);
dy_dt = f-exch;
end