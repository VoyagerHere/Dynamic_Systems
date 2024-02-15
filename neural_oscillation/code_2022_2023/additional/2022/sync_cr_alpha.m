sigma = linspace(0, 0.5);



g1 = 1.01;

alpha = pi;
delta1 = 0.15;
delta2 = 0.25;
delta3 = 1;
d1 = (abs(1)/2 + abs(delta1)./(4*sin(sigma/4)));
d2 = (abs(1)/2 + abs(delta2)./(4*sin(sigma/4)));
d3 = (abs(1)/2 + abs(delta3)./(4*sin(sigma/4)));
plot(sigma, d1, sigma, d2, sigma, d3);
ylim([0, 10]);
xlabel('\sigma'); 
ylabel('d_{cr}');
legend('|\sigma| = 0.15', '|\sigma| = 0.25', '|\sigma| = 1')
    