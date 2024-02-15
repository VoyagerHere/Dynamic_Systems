% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  Visualisator of phi-alpha system 
%
% WORK IN PROGRESS
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



d1 = GAMMA_TABLE(:, 1);
d2 = GAMMA_TABLE(:, 2);
ratio = GAMMA_TABLE(:, 3);
delta = GAMMA_TABLE(:, 4);

% Find unique ratios
unique_ratios = unique(ratio);

% Define colors for different ratios
colors = lines(length(unique_ratios));
gray_color = [0.8 0.8 0.8];


ratio(delta > 0 & ratio == 1) = -9;


% Plot points with same ratio in same color
figure;
hold on;

for i = 1:length(unique_ratios)
    ratio_index = find(ratio == unique_ratios(i));
    if unique_ratios(i) == -9
        % Plot points with ratio -5 in gray color
        plot3(d1(ratio_index), d2(ratio_index), delta(ratio_index), '.', 'Color', gray_color);
        legend_str = 'Quasi-Periodic';
    elseif unique_ratios(i) == -1 
        % Plot points with ratio -1 or -2 in darker color than gray
        plot3(d1(ratio_index), d2(ratio_index), delta(ratio_index), '.', 'Color', 0.4 * gray_color);
        legend_str = 'Death \phi_1';
    elseif unique_ratios(i) == -2
        plot3(d1(ratio_index), d2(ratio_index), delta(ratio_index), '.', 'Color', 0.6 * gray_color);
        legend_str = 'Death \phi_2';
    elseif unique_ratios(i) == -3
        plot3(d1(ratio_index), d2(ratio_index), delta(ratio_index), '.', 'Color', 'black');
        legend_str = 'Death \phi_3';
    elseif unique_ratios(i) == -4
        plot3(d1(ratio_index), d2(ratio_index), delta(ratio_index), '.', 'Color', 'black');
        legend_str = 'Death \phi_1, \phi_2, \phi_3';
    else
        % Plot points with other ratios in random colors
        plot3(d1(ratio_index), d2(ratio_index), delta(ratio_index), '.', 'Color', rand(1,3));
        legend_str = ['1:' num2str(unique_ratios(i))];
    end
    legend_entries{i} = legend_str;
end

lgd = legend(legend_entries);
fontsize(lgd,14,'points')
xlabel("d",'FontSize', 14,'FontWeight','bold');
ylabel('\Delta', 'Interpreter','tex','FontSize', 14,'FontWeight','bold');
title(sprintf('Δ max = %g, α = %s, n_1 = 2, n_2 = 2', delta_max, alpha_text));


% xlim([0 max(d)])
% ylim([0, max(delta)])