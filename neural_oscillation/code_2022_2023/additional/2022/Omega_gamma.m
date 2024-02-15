d = 0;
warning('off');

global n;
global g;
W1 = [];


gg = linspace(1.01, 2.01, 10);

for i = 1:length(gg)
    g = gg(i);
    disp(g);
    disp(n);
    w1 = f1();
    W1(i, 1) = w1;
end

plot(gg, W1(:,1));
xlabel('\gamma', 'Interpreter','tex');
ylabel('\Omega', 'Interpreter','tex');
xlim([min(gg), max(gg)]);

function w1 = f1()
  global g;
  w1 = sqrt(g^2-1);
end
