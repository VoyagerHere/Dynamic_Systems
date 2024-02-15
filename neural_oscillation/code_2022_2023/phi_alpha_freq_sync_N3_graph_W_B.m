freq_ratio = DiagOut(:, 1);
d = DiagOut(:, 3);
err = DiagOut(:, 5);


% freq_ratio = round(freq_ratio, 3);

freq_ratio = round(freq_ratio, 4);

unique_mode = unique(err);

gray_color = [0.8 0.8 0.8];


alpha_text = '2π/3';
delta_max = 0.01;

figure;
hold on;

for i = 1:length(unique_mode)
    mode_index = find(err == unique_mode(i));
    if unique_mode(i) == 1
        death_phi_1 = d(err == 1);
        area(death_phi_1, ones(size(death_phi_1)),'FaceColor', 0.6 * gray_color);
        legend_str = 'Death \phi_1';
    elseif unique_mode(i) == 2 
        death_phi_2 = d(err == 2);
        area(death_phi_2, ones(size(death_phi_2)),'FaceColor', 0.4 * gray_color);
        legend_str = 'Death \phi_2';
    elseif unique_mode(i) == 3
        death_phi_3 = d(err == 3);
        area(death_phi_3, ones(size(death_phi_3)),'FaceColor','black');
        legend_str = 'Death \phi_1, \phi_2';
    elseif unique_mode(i) == 4
        death_phi_4 = d(err == 4);
        area(death_phi_4, ones(size(death_phi_4)),'FaceColor','black');
        legend_str = 'Death \phi_1, \phi_2, \phi_3';
    else
        plot(d(err == 0), freq_ratio(err == 0), '.', 'Color', 'blue');
        legend_str = 'w_{b}';
    end
    legend_entries{i} = legend_str;
end
legend(legend_entries);
xlabel("d",'FontSize', 13,'FontWeight','bold');
ylabel('\omega_{b}^1/\omega_{b}^2', 'Interpreter','tex','FontSize', 13,'FontWeight','bold');
title(sprintf('Δ = %g, α = %s, n_1 = 3, n_2 = 3, n_3 = 3', delta_max, alpha_text));

hold off;