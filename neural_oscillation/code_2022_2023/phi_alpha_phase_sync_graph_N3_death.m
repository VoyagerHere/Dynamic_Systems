    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  Visualisation of phi-alpha system phase sync diagram
%  N = 2
%
% WORK IN PROGRESS
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


d = GAMMA_TABLE(:, 1);
d2 = GAMMA_TABLE(:, 2);
ratio = GAMMA_TABLE(:, 3);
ratio2 = GAMMA_TABLE(:, 4);
delta = GAMMA_TABLE(:, 5);
sync = GAMMA_TABLE(:, 6);
err = GAMMA_TABLE(:, 7);


% Find unique ratios
unique_ratios = unique(ratio);

% Define colors for different ratios
colors = lines(length(unique_ratios));
gray_color = [0.8 0.8 0.8];

% Plot points with same ratio in same color
figure;
hold on;

% blue gradient
left_color = [1 1 1]; % white
right_color = [0 0.25 0.5]; % color according to the attached image
cmap_blue = interp1([0, 1], [left_color; right_color], linspace(0, 1, 4));


legend_entries = {};
base_color = [0.3010 0.7450 0.9330];
blue_color = [0 0.4470 0.7410];

death_all = find(err == 111);
death_1_2 = find(err == 110);
death_1 = find(err == 100);
death_none = find(err == 0);
if(length(death_none) > 1)
    plot(d(death_none), delta(death_none), '.', 'Color', gray_color);
    legend_str = 'Quasi-Periodic';
    legend_entries{end+1} = legend_str;
end
if(length(death_1) > 1)
    plot(d(death_1), delta(death_1), '.', 'Color', 0.4 * gray_color);
    legend_str = 'Death \phi_1';
    legend_entries{end+1} = legend_str;
end
if(length(death_1_2) > 1)
    % Plot points with ratio -1 or -2 in darker color than gray
    plot(d(death_1_2), delta(death_1_2), '.', 'Color', 0.6 * gray_color);
    legend_str = 'Death \phi_1, \phi_2';
    legend_entries{end+1} = legend_str;
end
if(length(death_all) > 1)
    % Plot points with ratio -1 or -2 in darker color than gray
    plot(d(death_all), delta(death_all), '.', 'Color', 'black');
    legend_str = 'Death \phi_1, \phi_2, \phi_3';
    legend_entries{end+1} = legend_str;
end

lgd = legend(legend_entries);
fontsize(lgd,14,'points')
xlabel("d",'FontSize', 14,'FontWeight','bold');
ylabel('\Delta', 'Interpreter','tex','FontSize', 14,'FontWeight','bold');
title(sprintf('Δ max = %g, α = %s, n_1 = 1, n_2 = 1, n_3 = 1', delta_max, alpha_text));