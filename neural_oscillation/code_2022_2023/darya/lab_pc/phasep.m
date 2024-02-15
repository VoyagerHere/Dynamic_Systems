d=0.15;
sigma = 1;

n1 = 1;
n2 = 1;

y0 = [(pi/2); (pi/2)]; %временные рамки для построения
%y0 = [pi; pi];
a=10000;
b=12000;

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
N = min([length(T1),length(T2),length(T3),length(T4)]);
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

yFinal1 = [];
yFinal2 = [];
cutCount1 = 1;
cutStart1 = 1;
cutEnd1 = 1;

cutCount2 = 1;
cutStart2 = 1;
cutEnd2 = 1;

cutPoint = -1;

hold on
for i = 1:length(Y(:,1))
    if(abs(Y(i+1,2) - Y(i,2)) > (2*pi - 0.5))
        cutEnd2 = i; 
        yFinal2 = Y(cutStart2:cutEnd2,2);
        for j = cutStart2:cutEnd2
            if(abs(Y(j+1,1) - Y(j,1)) > (2*pi - 0.5))
                cutPoint = j;
                break
            end
        end
        if(cutPoint ~= -1)
            plot(Y(cutStart2:cutPoint - 1,1),yFinal2(1 : cutPoint - cutStart2),'-r');
            plot(Y(cutPoint + 2:cutEnd2,1),yFinal2(cutPoint - cutStart2 + 3 : end),'-r');
        else
            plot(Y(cutStart2:cutEnd2,1),yFinal2,'-r');
        end
        cutPoint = -1;
        cutStart2 = cutEnd2 + 2;
    end
end
hold off
%построение графика
%plot(Y(:,1), Y(:,2),'-r');

ylim([0 2*pi]); yticks([0 pi 2*pi]); yticklabels(["0" "\pi" "2\pi"]);
xlim([0 2*pi]); xticks([0 pi 2*pi]); xticklabels(["0" "\pi" "2\pi"]);
title('d=0.0015', 'sigma = 1');
legend('n_1=1','n_2=1', 'Interpreter','tex');
xlabel('\theta_1', 'Interpreter','tex');
ylabel('\theta_2', 'Interpreter','tex');
%xlabel('\phi_1', 'Interpreter','tex');
%ylabel('\phi_2', 'Interpreter','tex');

%задание ДУ

function dy_dt = eqn(~,y,F,d,n1,n2)
g = [1.01; 1.02];
no = [n1; n2];
f = g - sin(y./no);
exch = d*F([2;1]);
dy_dt = f-exch;
end
