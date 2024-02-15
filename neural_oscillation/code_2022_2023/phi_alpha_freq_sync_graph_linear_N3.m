freq_ratio = DiagOut(:, 2);
d = DiagOut(:, 3);
err = DiagOut(:, 5);



unique_mode = unique(err);

gray_color = [0.8 0.8 0.8];


alpha_text = 'π/8';
delta_max = 0.01;

figure;
hold on;


freq_ratio = round(freq_ratio, 2);
y_values = unique(freq_ratio); 



for i = 1:length(unique_mode)
    mode_index = find(err == unique_mode(i));
    if unique_mode(i) == 1
        death_phi_1 = d(err == 1);
        area(death_phi_1, ones(size(death_phi_1)),'FaceColor', 0.4 * gray_color);
        legend_str = 'Death \phi_1';
    elseif unique_mode(i) == 2 
        death_phi_2 = d(err == 2);
        area(death_phi_2, ones(size(death_phi_2)),'FaceColor', 0.6 * gray_color);
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
        for j = 1:length(y_values)
            value = y_values(j);
            if (isnan(value)) 
              continue;
            end
            x_values = d(freq_ratio == value);
%             if(length(x_values) < 4)
%                 continue;
%             end
            plot(x_values, value*ones(size(x_values)), 'b', 'LineWidth', 2);
        end
        legend_str = 'w_{s}';
    end
    legend_entries{i} = legend_str;
end
legend(legend_entries);
xlabel("d",'FontSize', 13,'FontWeight','bold');
ylabel('\omega_1/\omega_2', 'Interpreter','tex','FontSize', 13,'FontWeight','bold');
title(sprintf('Δ = %g, α = %s, n_1 = 3, n_2 = 3, n_3 = 3', delta_max, alpha_text));

hold off;