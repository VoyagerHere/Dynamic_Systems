x = linspace(-20,20,400);
n = 10;
y = sin(x/n);
plot(x,y)
grid on; grid minor;
xticks([-6*pi -5*pi -4*pi -3*pi -2*pi -pi 0 pi 2*pi 3*pi 4*pi 5*pi 6*pi])
xticklabels({'-6\pi','-5\pi', '-4\pi',  '-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi','4\pi', '5\pi', '6\pi'})
title('F(\phi), n=2', 'Interpreter','tex');
ylabel("y");
xlabel("x");
