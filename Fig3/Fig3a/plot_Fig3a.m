addpath(genpath('.\ccode'));
clear; clc; close all;

%% Plot the background phase diagram
chi11 = -7.0; chi22 = -7.0; chi12 = -4.5;
figure('Position', [150, 150, 800, 600], 'Color', 'w');
plot_cartesian_phase_diagram(chi11, chi12, chi22, [], 0.05);
hold on;

%% Plot the instability boundary
load('instability_boundary.mat')
critical_x = critical_points(:,1);
critical_y = critical_points(:,2);
plot(critical_x, critical_y, '-o', 'LineWidth', 2, 'MarkerSize', 4, ...
    'MarkerFaceColor','b', 'Color', 'b')

%% Plot simulation points
load('phase_points.mat')
color_homo  = [31, 119, 180] / 255;  
color_two   = [255, 127, 14] / 255;  
color_three = [214, 39, 40] / 255;  

colors  = {color_homo, color_two, color_three};
markers = {'o', 's', '^'}; 
labels  = {'Homogeneous', 'Two-phase', 'Three-phase'};

% Group points by the number of phases
for phase_num = 1:3
    idx = (phase == phase_num);
    if any(idx)
        scatter(phi1(idx), phi2(idx), 150, markers{phase_num}, ...
            'MarkerFaceColor', colors{phase_num}, ...
            'MarkerEdgeColor', 'k', ...   
            'LineWidth', 0.8, ...
            'DisplayName', labels{phase_num});
    end
end

xlim([0,0.32])
ylim([0,0.22])

set(gca, 'Layer', 'top', 'FontSize', 20, 'LineWidth', 1.2, 'TickDir', 'in');

hx = xlabel('$\bar\phi_1$', 'Interpreter', 'latex', 'FontSize', 25);
hy = ylabel('$\bar\phi_2$', 'Interpreter', 'latex', 'FontSize', 25);
ht = title('$\chi_{12}=-4.5$', 'Interpreter', 'latex', 'FontSize', 25);

set(ht, 'Units', 'normalized', 'Position', [0.5, 1.04, 0]);  

x_marks = [0.165669, 0.147694, 0.108, 0.120369];
y_marks = [0.0721281, 0.0951886, 0.16658, 0.20321];

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

    quiver(x_text, y_text, -dx, -dy, 0, ...
        'Color', mark_colors(i,:), ...
        'LineWidth', 2, ...   
        'MaxHeadSize', 4);     

    text(x_text, y_text, sprintf('(%d)', i), ...
         'FontSize', 16, ...
         'FontWeight', 'bold', ...
         'Color', mark_colors(i,:), ...
         'HorizontalAlignment', 'center', ... 
         'VerticalAlignment', 'middle', ...   
         'BackgroundColor', [1 1 1 0.85], ... 
         'Margin', 1);
end

set(gca, 'TickLabelInterpreter', 'latex', 'GridColor', [0.8 0.8 0.8]);

%% Legend and export
h_three = plot(NaN, NaN, '^', 'MarkerSize', 14, 'MarkerFaceColor', color_three, ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'LineStyle', 'none'); 
h_two   = plot(NaN, NaN, 's', 'MarkerSize', 14, 'MarkerFaceColor', color_two, ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'LineStyle', 'none');

lgd = legend([h_three, h_two], {'Three-phase', 'Two-phase'}, ...
    'Location', 'west', ...    
    'FontSize', 18, ...        
    'Box', 'on');
set(lgd, 'Interpreter', 'latex');

% exportgraphics(gcf, 'Add_suddenly.pdf', 'ContentType', 'vector');