sigma = SIGMA_TABLE(:, 1);
d = SIGMA_TABLE(:, 2); % if d==0, may be it non-sync process
ratio = SIGMA_TABLE(:, 3);
% phase_mode = GAMMA_TABLE(:, 4);
% 
% indices = find(phase_mode == -1);
% ratio(indices) = -3;
% indices_2 = find(ratio == -2);
% ratio(indices_2) = -1;

% Find unique ratios
unique_ratios = unique(ratio);

% Define colors for different ratios
colors = lines(length(unique_ratios));

% Plot points with same ratio in same color
figure;
hold on;

for i = 1:length(unique_ratios)
  ratio_index = find(ratio == unique_ratios(i));
  color_index = find(unique_ratios(i) == [-9, -2, -1, -3]) + 1;
  if unique_ratios(i) == -9
      legend_str = 'Quasi-Periodic';
      plot(sigma(ratio_index), d(ratio_index), '.', 'Color', [0 0.5 0]);
  elseif unique_ratios(i) == -1
      legend_str = 'Death \phi_1';
      plot(sigma(ratio_index), d(ratio_index), '.', 'Color', 	[0 0.4470 0.7410]);
  elseif unique_ratios(i) == -2
      legend_str = 'Death \phi_2';
      plot(sigma(ratio_index), d(ratio_index), '.', 'Color', 	[0.3010 0.7450 0.9330]);
  elseif unique_ratios(i) == -3
      legend_str = 'No oscillation';
      plot(sigma(ratio_index), d(ratio_index), '.', 'Color', 	[0, 0, 1]);
  else
      legend_str = ['1:' num2str(unique_ratios(i))];
      plot(sigma(ratio_index), d(ratio_index), '.', 'Color', rand(1,3));
  end
  legend_entries{i} = legend_str;
end

xlabel('σ');
ylabel('d');
title('Points with Same Ratio in Same Color');
legend(legend_entries);





% xlim([0 max(d)])
% ylim([0, max(delta)])