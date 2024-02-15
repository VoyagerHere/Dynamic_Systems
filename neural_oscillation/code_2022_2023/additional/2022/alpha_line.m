% create figure window
figure
hold on

% draw horizontal line
line([0, 2*pi], [0, 0], 'LineWidth', 1);

% draw dots at 0 and 2pi
plot([0, 2*pi], [0, 0], 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'k');

% draw vertical line at pi/2
line([pi/2, pi/2], [-0.8, 0.8], 'LineStyle', '--', 'LineWidth', 2);

% draw labels
text(0, -0.15, '0');
text(pi/2, -0.8, '\pi/2');
text(2*pi, -0.15, '2\pi');

% draw text
text(0.75, -1, 'Отталкивающая связь');
text(4.0, -1, 'Притягивающая связь');

% set axis limits and labels
xlim([-0.5, 2*pi+0.5]);
ylim([-1.2, 1.2]);
xlabel('Угол');
ylabel('Связь');

% turn off box and grid
box off
grid off

% hold off
hold off
