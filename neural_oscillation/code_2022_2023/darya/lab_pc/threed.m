d1=0.0015;
d2=0.0015;
sigma = 1;

n1 = 1;
n2 = 1;
n3 = 1;

y0 = [(2*pi); (2*pi);(2*pi)]; %временные рамки для построения

a=10000;
b=12000;

opts = odeset('RelTol',2e-13,'AbsTol',1e-100); %решение системы ДУ
F = [1 ; 1; 1];
[T1,Y1] = ode45(@(t,y)eqn(t,y,F,d1, d2,n1,n2, n3), [a, b], y0, opts);

F = [0 ; 1; 0];
[T2,Y2] = ode45(@(t,y)eqn(t,y,F,d1, d2,n1,n2, n3), [a, b], y0, opts);

F = [0 ; 1; 1];
[T3,Y3] = ode45(@(t,y)eqn(t,y,F,d1, d2,n1,n2, n3), [a, b], y0, opts);

F = [0 ; 0; 1];
[T4,Y4] = ode45(@(t,y)eqn(t,y,F,d1, d2,n1,n2, n3), [a, b], y0, opts);

F = [0 ; 0; 0];
[T5,Y5] = ode45(@(t,y)eqn(t,y,F,d1, d2,n1,n2, n3), [a, b], y0, opts);

F = [1 ; 0; 0];
[T6,Y6] = ode45(@(t,y)eqn(t,y,F,d1, d2,n1,n2, n3), [a, b], y0, opts);

F = [1 ; 0; 1];
[T7,Y7] = ode45(@(t,y)eqn(t,y,F,d1, d2,n1,n2, n3), [a, b], y0, opts);

F = [1 ; 1; 0];
[T8,Y8] = ode45(@(t,y)eqn(t,y,F,d1, d2,n1,n2, n3), [a, b], y0, opts);




%блок проверки решений на соответствие принадлежности к D_i

for(i = 1:length(T1)) 
    if((Y1(i,1) < (3*pi/2) + sigma && Y1(i,1) > (3*pi/2) - sigma) || (Y1(i,2) < (3*pi/2) + sigma && Y1(i,2) > (3*pi/2) - sigma)|| (Y1(i,3) < (3*pi/2) + sigma && Y1(i,3) > (3*pi/2) - sigma))
        Y1(i,1) = -1;
        Y1(i,2) = -1;
        Y1(i,3) = -1;
    end
end

for(i = 1:length(T2))
    if((Y2(i,1) > (3*pi/2) + sigma || Y2(i,1) < (3*pi/2) - sigma) || (Y2(i,2) < (3*pi/2) + sigma && Y2(i,2) > (3*pi/2) - sigma) ||(Y2(i,3) > (3*pi/2) + sigma || Y2(i,3) < (3*pi/2) - sigma))
        Y2(i,1) = -1;
        Y2(i,2) = -1;
        Y2(i,3) = -1;
    end
end

for(i = 1:length(T3))
    if((Y3(i,1) > (3*pi/2) + sigma || Y3(i,1) < (3*pi/2) - sigma) || (Y3(i,2) < (3*pi/2) + sigma && Y3(i,2) > (3*pi/2) - sigma) || (Y3(i,3) < (3*pi/2) + sigma && Y3(i,3) > (3*pi/2) - sigma))
        Y3(i,1) = -1;
        Y3(i,2) = -1;
        Y3(i,3) = -1;
    end
end

for(i = 1:length(T4))
    if((Y4(i,1) > (3*pi/2) + sigma || Y4(i,1) < (3*pi/2) - sigma) || (Y4(i,2) > (3*pi/2) + sigma || Y4(i,2) < (3*pi/2) - sigma) ||(Y4(i,3) < (3*pi/2) + sigma && Y4(i,3) > (3*pi/2) - sigma))
        Y4(i,1) = -1;
        Y4(i,2) = -1;
        Y4(i,3) = -1;
    end
end

for(i = 1:length(T5))
    if((Y5(i,1) > (3*pi/2) + sigma || Y5(i,1) < (3*pi/2) - sigma) || (Y5(i,2) > (3*pi/2) + sigma || Y5(i,2) < (3*pi/2) - sigma) ||(Y5(i,3) > (3*pi/2) + sigma || Y5(i,3) < (3*pi/2) - sigma))
        Y5(i,1) = -1;
        Y5(i,2) = -1;
        Y5(i,3) = -1;
    end
end

for(i = 1:length(T6))
    if((Y6(i,1) < (3*pi/2) + sigma && Y6(i,1) > (3*pi/2) - sigma) || (Y6(i,2) > (3*pi/2) + sigma || Y6(i,2) < (3*pi/2) - sigma) ||(Y6(i,3) > (3*pi/2) + sigma || Y6(i,3) < (3*pi/2) - sigma))
        Y6(i,1) = -1;
        Y6(i,2) = -1;
        Y6(i,3) = -1;
    end
end

for(i = 1:length(T7))
    if((Y7(i,1) < (3*pi/2) + sigma && Y7(i,1) > (3*pi/2) - sigma) || (Y7(i,2) > (3*pi/2) + sigma || Y7(i,2) < (3*pi/2) - sigma) ||(Y7(i,3) < (3*pi/2) + sigma && Y7(i,3) > (3*pi/2) - sigma))
        Y7(i,1) = -1;
        Y7(i,2) = -1;
        Y7(i,3) = -1;
    end
end

for(i = 1:length(T8)) 
    if((Y8(i,1) < (3*pi/2) + sigma && Y8(i,1) > (3*pi/2) - sigma) || (Y8(i,2) < (3*pi/2) + sigma && Y8(i,2) > (3*pi/2) - sigma)|| (Y8(i,3) > (3*pi/2) + sigma || Y8(i,3) < (3*pi/2) - sigma))
        Y8(i,1) = -1;
        Y8(i,2) = -1;
        Y8(i,3) = -1;
    end
end

%внесение всех подходящих решений в общий массив

Y = [];
N = min([length(T1),length(T2),length(T3),length(T4),length(T5),length(T6),length(T7),length(T8)]);
for(i = 1:N)
    if(Y1(i,1) > -1)
        Y(i,1) = Y1(i,1);
        Y(i,2) = Y1(i,2);
        Y(i,3) = Y1(i,3);
    end
    if(Y2(i,1) > -1)
        Y(i,1) = Y2(i,1);
        Y(i,2) = Y2(i,2);
        Y(i,3) = Y2(i,3);
    end
    if(Y3(i,1) > -1)
        Y(i,1) = Y3(i,1);
        Y(i,2) = Y3(i,2);
        Y(i,3) = Y3(i,3);
    end
    if(Y4(i,1) > -1)
        Y(i,1) = Y4(i,1);
        Y(i,2) = Y4(i,2);
        Y(i,3) = Y4(i,3);
    end
    if(Y5(i,1) > -1)
        Y(i,1) = Y5(i,1);
        Y(i,2) = Y5(i,2);
        Y(i,3) = Y5(i,3);
    end
    if(Y6(i,1) > -1)
        Y(i,1) = Y6(i,1);
        Y(i,2) = Y6(i,2);
        Y(i,3) = Y6(i,3);
    end
    if(Y7(i,1) > -1)
        Y(i,1) = Y7(i,1);
        Y(i,2) = Y7(i,2);
        Y(i,3) = Y7(i,3);
    end
    if(Y8(i,1) > -1)
        Y(i,1) = Y8(i,1);
        Y(i,2) = Y8(i,2);
        Y(i,3) = Y8(i,3);
    end
end

Y = mod(Y, 2*pi);

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
if(N == length(T5))
    T = T5;
end
if(N == length(T6))
    T = T6;
end
if(N == length(T7))
    T = T7;
end
if(N == length(T8))
    T = T8;
end

yFinal1 = [];
yFinal2 = [];
yFinal3 = [];
cutCount1 = 1;
cutStart1 = 1;
cutEnd1 = 1;

cutCount2 = 1;
cutStart2 = 1;
cutEnd2 = 1;

cutCount3 = 1;
cutStart3 = 1;
cutEnd3 = 1;

cutPoint1 = -1;
cutPoint2 = -1;

hold on
for i = 1:length(Y(:,1))
    if(i+1 < length(Y(:,1)))
        if(abs(Y(i+1,3) - Y(i,3)) > (2*pi - 0.5))
            cutEnd3 = i; 
            yFinal3 = Y(cutStart3:cutEnd3,3);
            for j = cutStart3:cutEnd3
                if(abs(Y(j+1,1) - Y(j,1)) > (2*pi - 0.5))
                    cutPoint1 = j;
                    break
                end
                if(abs(Y(j+1,2) - Y(j,2)) > (2*pi - 0.5))
                    cutPoint2 = j;
                    break
                end
            end
            if(cutPoint1 ~= -1 && cutPoint2 ~= -1)
                if(cutPoint1 > cutPoint2)
                    plot3(Y(cutStart3:cutPoint1 - 1,1),Y(cutStart3:cutPoint1 - 1,2),yFinal3(1 : cutPoint1 - cutStart3),'-r');
                    plot3(Y(cutPoint1 + 2:cutPoint2 - 1,1),Y(cutPoint1 + 2:cutPoint2 - 1,2),yFinal3(cutPoint1 - cutStart3 + 1 : cutPoint2 - cutStart3),'-r');
                    plot3(Y(cutPoint2 + 2:cutEnd3,1),Y(cutPoint2 + 2:cutEnd3,2),yFinal3(cutPoint2 - cutStart3 + 3: end),'-r');
                else
                    plot3(Y(cutStart3:cutPoint2 - 1,1),Y(cutStart3:cutPoint2 - 1,2),yFinal3(1 : cutPoint2 - cutStart3),'-r');
                    plot3(Y(cutPoint2 + 2:cutPoint1 - 1,1),Y(cutPoint2 + 2:cutPoint1 - 1,2),yFinal3(cutPoint2 - cutStart3 + 1 : cutPoint1 - cutStart3),'-r');
                    plot3(Y(cutPoint1 + 2:cutEnd3,1),Y(cutPoint1 + 2:cutEnd3,2),yFinal3(cutPoint1 - cutStart3 + 3: end),'-r');
                end
            elseif(cutPoint1 ~= -1 && cutPoint2 == -1)
                plot3(Y(cutStart3:cutPoint1 - 1,1),Y(cutStart3:cutPoint1 - 1,2),yFinal3(1 : cutPoint1 - cutStart3),'-r');
                plot3(Y(cutPoint1 + 2:cutEnd3,1),Y(cutPoint1 + 2:cutEnd3,2),yFinal3(cutPoint1 - cutStart3 + 3 : end),'-r');
            elseif(cutPoint1 == -1 && cutPoint2 ~= -1)
                plot3(Y(cutStart3:cutPoint2 - 1,1),Y(cutStart3:cutPoint2 - 1,2),yFinal3(1 : cutPoint2 - cutStart3),'-r');
                plot3(Y(cutPoint2 + 2:cutEnd3,1),Y(cutPoint2 + 2:cutEnd3,2),yFinal3(cutPoint2 - cutStart3 + 3 : end),'-r');
            elseif(cutPoint1 == -1 && cutPoint2 == -1)
                plot3(Y(cutStart3:cutEnd3,1),Y(cutStart3:cutEnd3,2),yFinal3,'-r');
            end
            cutPoint1 = -1;
            cutPoint2 = -1;
            cutStart3 = cutEnd3 + 2;
        end
    end
end
hold off

view(3);

%построение графика ФП

%plot3(Y(:,1), Y(:,2), Y(:,3),'-r');
ylim([0 2*pi]); yticks([0 pi 2*pi]); yticklabels(["0" "\pi" "2\pi"]);
xlim([0 2*pi]); xticks([0 pi 2*pi]); xticklabels(["0" "\pi" "2\pi"]);
zlim([0 2*pi]); zticks([0 pi 2*pi]); zticklabels(["0" "\pi" "2\pi"]);
title('d_1=d_2=0.0015', 'sigma = 1');
legend('n_1=1','n_2=1','n_3=1', 'Interpreter','tex');
xlabel('\theta_1', 'Interpreter','tex');
ylabel('\theta_2', 'Interpreter','tex');
zlabel('\theta_3', 'Interpreter','tex');

%построение графика временной реализации
%plot(T, Y(:,1),'-b', T, Y(:,2),'-r', T, Y(:,3),'-g');
%ylim([0 2*pi]); yticks([0 pi 2*pi]); yticklabels(["0" "\pi" "2\pi"]);
%title('d = 0.0015', 'sigma = 1');
%legend('n_1=1','n_2=1', 'n_3=1', 'Interpreter','tex');
%xlabel("t");
%ylabel('\theta', 'Interpreter','tex');

%задание ДУ

function dy_dt = eqn(~,y,F,d1,d2,n1,n2,n3)
g = [1.01; 1.02;1.03];
no = [n1; n2; n3];
f = g - cos(y./no);
exch = d1*F([3;1;2])+d2*F([2;3;1]);
dy_dt = f-exch;
end
