addpath(genpath('.\ccode'));
clear; clc; close all;

%% Plot background phase diagram
chi11 = -7.0; chi22 = -7.0; chi12 = -4.5;
figure('Color', 'w', 'Position', [200, 100, 800, 550]);
plot_cartesian_phase_diagram(chi11, chi12, chi22, [], 0.05);
hold on;

%% Plot instability boundary
load('instability_boundary.mat')
critical_x = critical_points(:,1);
critical_y = critical_points(:,2);
plot(critical_x, critical_y, '-o', 'Color', 'b', 'MarkerSize', 4, ...
    'MarkerFaceColor', 'b', 'LineWidth', 2);
plot(critical_y, critical_x, '-o', 'Color', 'b', 'MarkerSize', 4, ...
    'MarkerFaceColor','b', 'LineWidth', 2);

%% Plot simulation phase points
load('phase_points2.mat')
color_homo  = [31, 119, 180] / 255;  
color_two   = [255, 127, 14] / 255;  
color_three = [214, 39, 40] / 255;  
colors  = {color_homo, color_two, color_three};
markers = {'o', 's', '^'};
labels  = {'Homogeneous', 'Two-phase', 'Three-phase'};

% Plot different symbols according to the number of phases
for phase_num = 1:3
    idx = (phase == phase_num);
    if any(idx)
        scatter(phi1(idx), phi2(idx), 80, markers{phase_num}, ...
            'MarkerFaceColor', colors{phase_num}, ...
            'MarkerEdgeColor', 'k', ...   
            'LineWidth', 0.8);
    end
end

%% Representative state
h_dummy = gobjects(1, 3); 

for k = 1:3
    h_dummy(k) = plot(NaN, NaN, markers{k}, ...
        'MarkerFaceColor', colors{k}, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 10, ...       
        'LineWidth', 0.8, ...
        'LineStyle', 'none');      
end

lgd = legend(h_dummy, labels, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 16);
set(lgd, 'Box', 'on'); 

set(gca, 'Position', [0.15, 0.15, 0.75, 0.75]); 
set(gca, 'Layer', 'top', 'FontSize', 20, 'LineWidth', 1.2, 'TickDir', 'in');
box on

xlim([0.05, 0.26])
ylim([0.05, 0.26])
set(gca, 'XTick', 0.06:0.06:0.24);
set(gca, 'YTick', 0.06:0.06:0.24);

hy = ylabel(sprintf('$\\bar{\\phi}_2 $'), 'Interpreter', 'latex', 'FontSize', 25);
hx = xlabel(sprintf('$\\bar{\\phi}_1 $'), 'Interpreter', 'latex', 'FontSize', 25);
ht = title(sprintf('$\\chi_{12}=-4.5 $'), 'Interpreter', 'latex', 'FontSize', 25);

set(hx, 'Units', 'normalized', 'Position', [0.5, -0.09, 0]); 
set(hy, 'Units', 'normalized', 'Position', [-0.15, 0.5, 0]);
set(ht, 'Units', 'normalized', 'Position', [0.5, 1.03, 0]);

x_marks = [0.105873, 0.160342, 0.179262, 0.107518];
y_marks = [0.095714, 0.077886, 0.111275, 0.165370];

mark_colors = [
    1.00, 0.00, 0.00;  
    0.00, 0.65, 0.00;  
    0.00, 0.45, 0.85;  
    0.60, 0.20, 0.80   
];

dx = 0.017; 
dy = 0.000;
for i = 1:length(x_marks)
    x_target = x_marks(i);
    y_target = y_marks(i);
    
    x_text = x_target + dx;
    y_text = y_target + dy;

    h_arrow =quiver(x_text, y_text, -dx, -dy, 0, ...
        'Color', mark_colors(i,:), ...
        'LineWidth', 2, ...       
        'MaxHeadSize', 4);       
    h_arrow.Annotation.LegendInformation.IconDisplayStyle = 'off';

    text(x_text, y_text, sprintf('(%d)', i), ...
         'FontSize', 16, ...
         'FontWeight', 'bold', ...
         'Color', mark_colors(i,:), ...
         'HorizontalAlignment', 'center', ... 
         'VerticalAlignment', 'middle', ...   
         'BackgroundColor', [1 1 1 0.85], ... 
         'Margin', 1);
end
box on;

%% export
% exportgraphics(gcf, 'Add_Simultaneously_chi12.pdf', 'ContentType', 'vector');