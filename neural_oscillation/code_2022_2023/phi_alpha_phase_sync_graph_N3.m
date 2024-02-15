% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  Visualisation of phi-alpha system phase sync diagram
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

d = GAMMA_TABLE(:, 1);
d2 = GAMMA_TABLE(:, 2);
ratio1 = GAMMA_TABLE(:, 3);
ratio2 = GAMMA_TABLE(:, 4);
delta = GAMMA_TABLE(:, 5);
sync = GAMMA_TABLE(:, 6);




% Find unique ratios
unique_ratios = unique(ratio1);

% Define colors for different ratios
colors = lines(length(unique_ratios));
gray_color = [0.8 0.8 0.8];

% Plot points with same ratio in same color
figure;
hold on;

ylim([0 max(delta)]);

% blue gradient
left_color = [1 1 1]; % white
right_color = [0 0.25 0.5]; % color according to the attached image
cmap_blue = interp1([0, 1], [left_color; right_color], linspace(0, 1, 4));

legend_entries = {};
base_color = [0.3010 0.7450 0.9330];
blue_color = [0 0.4470 0.7410];





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% NOISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sync_full_noise = find((sync == 111));
sync_full_noise = intersect(sync_full_noise, find(delta>0));
sync(sync_full_noise) = 100;


% MIX BIG RATIO
% ind_1 = find(ratio1 > 3);
% ind_2 = find(ratio2 > 3);
% ratio1(ind_1) = 10;
% ratio2(ind_2) = 10;


% DELETE RECORDS
% indices = unique_ratios < 1;
% unique_ratios(indices) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


QUASI_PERIODIC = find(sync == 100);
plot(d(QUASI_PERIODIC), delta(QUASI_PERIODIC), '.', 'Color', gray_color);
legend_str = 'Quasi-Periodic';
legend_entries{end+1} = legend_str;

DEATH_ALL = find(ratio1 == -111);
plot(d(DEATH_ALL), delta(DEATH_ALL), '.', 'Color', 'black');
if ~(isempty(DEATH_ALL))
    legend_str = 'Death \phi_1, \phi_2, \phi_3';
    legend_entries{end+1} = legend_str;
end



% 
% BURST_FULL = find(sync == 101);
% plot(d(BURST_FULL), delta(BURST_FULL), '.', 'Color', rand(1,3));
% if ~(isempty(BURST_FULL))
%     legend_str = 'Burst sync global';
%     legend_entries{end+1} = legend_str;
% end

BURST_1_2 = find(sync == 202);
plot(d(BURST_1_2), delta(BURST_1_2), '.', 'Color', rand(1,3));
if ~(isempty(BURST_1_2))
  legend_str = 'Burst sync \phi_1, \phi_2';
  legend_entries{end+1} = legend_str;
end

BURST_2_3 = find(sync == 201);
plot(d(BURST_2_3), delta(BURST_2_3), '.', 'Color', rand(1,3));
if ~(isempty(BURST_2_3))
  legend_str = 'Burst sync \phi_2, \phi_3';
  legend_entries{end+1} = legend_str;
end

SPIKE_FULL = find(sync == 111);
plot(d(SPIKE_FULL), delta(SPIKE_FULL), '.', 'Color', rand(1,3));
if ~(isempty(SPIKE_FULL))
    legend_str = 'Spike sync global';
    legend_entries{end+1} = legend_str;
end

SPIKE_1_2 = find(sync == 212);
plot(d(SPIKE_1_2), delta(SPIKE_1_2), '.', 'Color', rand(1,3));
if ~(isempty(SPIKE_1_2))
  legend_str = 'Spike sync \phi_1, \phi_2';
  legend_entries{end+1} = legend_str;
end

SPIKE_2_3 = find(sync == 211);
plot(d(SPIKE_2_3), delta(SPIKE_2_3), '.', 'Color', rand(1,3));
if ~(isempty(SPIKE_2_3))
  legend_str = 'Spike sync \phi_2, \phi_3';
  legend_entries{end+1} = legend_str;
end



lgd = legend(legend_entries);
% objhl = findobj(lgd, 'type', 'line'); %// objects of legend of type line
% set(objhl, 'Markersize', 25); %// set marker size as desired
fontsize(lgd,14,'points')
xlabel("d",'FontSize', 14,'FontWeight','bold');
ylabel('\Delta', 'Interpreter','tex','FontSize', 14,'FontWeight','bold');
title(sprintf('Δ_{max} = %g, α = %s, n_1 = 3, n_2 = 3, n_3 = 3', delta_max, alpha_text),'FontSize', 14,'FontWeight','bold');  