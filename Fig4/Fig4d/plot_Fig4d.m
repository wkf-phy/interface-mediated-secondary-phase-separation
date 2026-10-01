addpath(genpath('.\ccode'));
load('interfacial_enrichment_kinetics.mat')

%% Plot kinetics
t_scaled = time_axis / 1e4;         
fig2 = figure('Color', 'w', 'Position', [300, 200, 700, 500]); 

% Total perimeter
yyaxis left
p1 = plot(t_scaled, Total_Perimeter, '-', 'LineWidth', 3.0, 'Color', [0 0.31 0.53]); 
hy = ylabel('Total Perimeter, $P_{\mathrm{total}}$', 'Interpreter', 'latex', 'FontSize', 20);

set(hy, 'Units', 'normalized', 'Position', [-0.1, 0.5, 0]);
set(gca, 'ycolor', [0 0.31 0.53], 'FontSize', 20, 'LineWidth', 1.5, 'TickDir', 'in');
ylim_left = ylim;
ylim([0, ylim_left(2)*1.05]); 

% Total adsorbed amount
yyaxis right
p2 = plot(t_scaled, Total_Adsorbed, '-', 'LineWidth', 3.0, 'Color', [0.85 0.33 0.10]); 

hr = ylabel('$M_{\mathrm{total}}$', 'Interpreter', 'latex', 'FontSize', 20);

ax = gca;
ax.YAxis(2).Exponent = 0; 

ylim_right = ylim;
ylim([0, ylim_right(2)*1.5]);

% Figure layout
set(gca, 'Position', [0.18, 0.18, 0.64, 0.68]);
xlim([0, 15]);
ax.XAxis.Exponent = 0; 
hx = xlabel('Time, $t$ ($\times 10^{4}$)', 'Interpreter', 'latex', 'FontSize', 20, 'FontWeight', 'bold');
set(hx, 'Units', 'normalized', 'Position', [0.5, -0.15, 0]);

ht = title('Kinetics of Interfacial Enrichment', 'Interpreter', 'latex', 'FontSize', 20);
set(ht, 'Units', 'normalized', 'Position', [0.5, 1.06, 0]);

% Grid and axes
grid on;
ax.GridColor = [0.1 0.1 0.1];
ax.GridAlpha = 0.15; 
ax.GridLineStyle = ':';
box on; 

%% Export
% exportgraphics(gcf, 'interfacial_enrichment_kinetics.pdf', 'ContentType', 'vector');
