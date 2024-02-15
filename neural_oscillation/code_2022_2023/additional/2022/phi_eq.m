
t = linspace(0, 600, 10000);
n1 = 2;
n2 = 3;
g = 1.01;
c = 0;
y1 = 2*n1*atan((1/g)*((1+sqrt(g^2-1)*tan(c+(sqrt(g^2-1)*t)/(2*n1)))));
y2 = 2*n2*atan((1/g)*((1+sqrt(g^2-1)*tan(c+(sqrt(g^2-1)*t)/(2*n2)))));

plot(t, y1,'-b',t, y2,'-g');
legend('n=2', 'n=3');
ylim([0 2*pi]); yticks([0 pi 3*pi/2 2*pi]); yticklabels(["0" "\pi" "3\pi/2" "2\pi"]); 
grid on; grid minor; 
title('\gamma = 1.01')
xlabel("t");
ylabel('\phi', 'Interpreter','tex');