% Load critical concentration profiles
load("critical_profiles.mat")

% Scan list and profile colors
c1_scan_list = [0.06, 0.10, 0.15, 0.20]; 
color_map = lines(length(c1_scan_list));

fig = figure('Position', [100, 100, 670, 500], 'Color', 'w'); 

%% Main panel: client concentration profiles
ax_main = axes('Position', [0.12, 0.12, 0.8, 0.8]);
hold(ax_main, 'on'); box(ax_main, 'on');

for i = 1:length(c1_scan_list)
    plot(ax_main, critical_profiles(i).r, critical_profiles(i).p2, ...
        'LineWidth', 2, 'Color', color_map(i,:), ...
        'DisplayName', sprintf('$\\bar{\\phi}_{1} = %.2f$', critical_profiles(i).c1_base));
end

legend(ax_main, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 15);

set(ax_main, 'TickLabelInterpreter', 'latex','FontSize', 14, 'LineWidth', 1.2, 'GridAlpha', 0.15);
xlabel(ax_main, '$r$', 'FontSize', 20, 'Interpreter', 'latex');
ylabel(ax_main, '$\phi_2$', 'FontSize', 20, 'Interpreter', 'latex');

legend(ax_main, 'Location', 'southeast', 'FontSize', 12);

ylim_main = ylim(ax_main);
yticks(ax_main, linspace(0,0.45,4));
title('Critical profiles')

%% Inset: scaffold concentration profiles
ax_inset = axes('Position', [0.68, 0.62, 0.2, 0.25]); 
hold(ax_inset, 'on'); box(ax_inset, 'on'); 

for i = 1:length(c1_scan_list)
    plot(ax_inset, critical_profiles(i).r, critical_profiles(i).p1, ...
        'LineWidth', 1.5, 'Color', color_map(i,:)); 
end

set(ax_inset, 'TickLabelInterpreter', 'latex','FontSize', 10, 'LineWidth', 1, 'GridAlpha', 0.15);
xlabel(ax_inset, '$r$', 'FontSize', 12, 'Interpreter', 'latex');
ylabel(ax_inset, '$\phi_1$', 'FontSize', 12, 'Interpreter', 'latex');

ylim_inset = ylim(ax_inset);
yticks(ax_inset, linspace(0,0.8, 3));

% exportgraphics(gcf, 'crtical_profiles.pdf', 'ContentType', 'vector');