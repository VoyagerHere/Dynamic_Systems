% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  Visualisation of phi-alpha system phase sync diagram
%  N = 2
%
% WORK IN PROGRESS
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



d = GAMMA_TABLE(:, 1);
ratio = GAMMA_TABLE(:, 2);
delta = GAMMA_TABLE(:, 3);
sync = GAMMA_TABLE(:, 4);

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



for i = 1:length(rmmissing(unique_ratios))
    ratio_index = find(ratio == unique_ratios(i));
    if unique_ratios(i) == 0
        % Plot points with ratio -5 in gray color
        plot(d(ratio_index), delta(ratio_index), '.', 'Color', gray_color);
        legend_str = 'S1: Quasi-Periodic';
        legend_entries{end+1} = legend_str;
    elseif unique_ratios(i) == -11 
        % Plot points with ratio -1 or -2 in darker color than gray
        plot(d(ratio_index), delta(ratio_index), '.', 'Color', 'black');
        legend_str = 'S4: Death \phi_1, \phi_2';
        legend_entries{end+1} = legend_str;
    elseif unique_ratios(i) == -10
        plot(d(ratio_index), delta(ratio_index), '.', 'Color', 0.6 * gray_color);
        legend_str = 'S3: Death \phi_2';
        legend_entries{end+1} = legend_str;
    elseif unique_ratios(i) == -1
        plot(d(ratio_index), delta(ratio_index), '.', 'Color', 0.4 * gray_color);
        legend_str = 'S2: Death \phi_1';
        legend_entries{end+1} = legend_str;
    else
        Burst = find(sync == 101);
        Spike = find(sync == 111);

        Burst_p = intersect(Burst, ratio_index);
        Spike_p = intersect(Spike, ratio_index);
        % Plot points with other ratios in random colors

        plot(d(Burst_p), delta(Burst_p), '.', 'Color', rand(1,3));
        plot(d(Spike_p), delta(Spike_p), '.', 'Color', rand(1,3));
        if ~(isempty(Burst_p))
            legend_str = ['Burst sync 1:' num2str(unique_ratios(i))];
            legend_entries{end+1} = legend_str;
        end
        if ~(isempty(Spike_p))
            legend_str = ['Spike sync 1:' num2str(unique_ratios(i))];
            legend_entries{end+1} = legend_str;
        end
    end
end

lgd = legend(legend_entries);
fontsize(lgd,14,'points')
xlabel("d",'FontSize', 14,'FontWeight','bold');
ylabel('\Delta', 'Interpreter','tex','FontSize', 14,'FontWeight','bold');
title(sprintf('Δ max = %g, α = %s, n_1 = 3, n_2 = 3', delta_max, alpha_text),'FontSize', 14,'FontWeight','bold');