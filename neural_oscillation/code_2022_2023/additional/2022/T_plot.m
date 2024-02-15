
% 
% LOW ACCURATE
% 
y0 = [ 0.0; 0.0];



% Y = mod(Y, 2*pi); 


% plot(T, Y(:,1),'-b');
% ylim([0 2*pi]); yticks([0 pi 2*pi 4*pi 6*pi 8*pi 10*pi]);
% yticklabels(["0" "\pi" "2\pi" "4\pi" "6\pi" "8\pi"  "10\pi"]);
OUT = [];
for i = 3:10
    n = i;
    opts = odeset('RelTol',2e-13,'AbsTol',1e-100, 'MaxStep',0.1);
    [T,Y] = ode45(@(t,y)eqn(t,y,n),[0, 6000], y0,opts);
    tilda = tld(n);

    
    for k = 1:n 
        mid1 = tilda+2*pi*(k-1);
        mid2 = mid1 + 2*pi;
        
        A1 = find(abs(Y(:,1) - mid1) < 0.04);
        A2 = find(abs(Y(:,1) - mid2) < 0.04);
        AA1 = mean(f4(A1,T));
        AA2 = mean(f4(A2,T));
        OUT(k, i-2) = AA2 - AA1;
    end
end
y3 = OUT(1:3,1);
y4 = OUT(1:4,2);
y5 = OUT(1:5,3);
y6 = OUT(1:6,4);
y7 = OUT(1:7,5);
y8 = OUT(1:8,6);
y9 = OUT(1:9,7);
y10 = OUT(1:10,8);


x3 = 1:1:3;
x4 = 1:1:4;
x5 = 1:1:5;
x6 = 1:1:6;
x7 = 1:1:7;
x8 = 1:1:8;
x9 = 1:1:9;
x10 = 1:1:10;



y3 = sort(y3,'descend');
y4 = sort(y4,'descend');
y5 = sort(y5,'descend');
y6 = sort(y6,'descend');
y7 = sort(y7,'descend');
y8 = sort(y8,'descend');
y9 = sort(y9,'descend');
y10 = sort(y10,'descend');


plot(x3,y3,x4,y4,x5,y5, x6, y6, x7, y7, x8, y8, x9, y9, x10, y10);

legend('n = 3', 'n = 4', 'n = 5', 'n = 6', 'n = 7', 'n = 8', 'n = 9', 'n = 10');
xlabel('T_n', 'Interpreter','tex');
ylabel("t");
title('\gamma = 1.01', 'Interpreter','tex');

ylim ([0 350]);

function AA1 = f4(AA1, T)
    AA1 = num2cell(AA1);
    AA1 = T([AA1{:}], :);
end

function tilda = tld(n)
    tilda = n*pi/2 - floor(n/4)*2*pi;
    if abs(tilda - pi/2) < 0.02
        tilda = 3*pi/2;
    else
        tilda = abs(tilda - pi);
    end
end


function dy_dt = eqn(t,y,n)
  g = 1.01;
  dy_dt = g-sin(y./n);
end
