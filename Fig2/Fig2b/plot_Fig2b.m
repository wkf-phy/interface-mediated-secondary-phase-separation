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
    'MarkerFaceColor', 'b', 'Color', 'b')

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
ylim([0,0.2])

set(gca, 'Layer', 'top', 'FontSize', 20, 'LineWidth', 1.2, 'TickDir', 'in', 'TickLabelInterpreter', 'latex');

hx = xlabel('$\phi_1$', 'Interpreter', 'latex', 'FontSize', 25);
hy = ylabel('$\phi_2$', 'Interpreter', 'latex', 'FontSize', 25);
ht = title('$\chi_{12}=-4.5$', 'Interpreter', 'latex', 'FontSize', 25);

set(ht, 'Units', 'normalized', 'Position', [0.5, 1.04, 0]);
box on;

%% Calculate the plot-box position for annotation arrows
drawnow; 

ax1 = gca;
fig = gcf;

fig.Units = 'pixels';
fig_W = fig.Position(3);
fig_H = fig.Position(4);

ax1.Units = 'normalized';
ax_pos = ax1.Position; 
x_lim = ax1.XLim;
y_lim = ax1.YLim;

pbar = ax1.PlotBoxAspectRatio;
box_ratio = pbar(1) / pbar(2);

avail_W_px = ax_pos(3) * fig_W;
avail_H_px = ax_pos(4) * fig_H;
avail_ratio = avail_W_px / avail_H_px; 

if box_ratio > avail_ratio
    true_W_norm = ax_pos(3);
    true_H_px = avail_W_px / box_ratio;
    true_H_norm = true_H_px / fig_H;
    true_X_norm = ax_pos(1);
    true_Y_norm = ax_pos(2) + (ax_pos(4) - true_H_norm) / 2;
else
    true_H_norm = ax_pos(4);
    true_W_px = avail_H_px * box_ratio;
    true_W_norm = true_W_px / fig_W;
    true_X_norm = ax_pos(1) + (ax_pos(3) - true_W_norm) / 2;
    true_Y_norm = ax_pos(2);
end

%% Mark selected parameter sets

x_marks = [ 0.15,0.15,];
y_marks = [ 0.10,0.08];

mark_colors = [
    0.60, 0.20, 0.80;  
    0.00, 0.80, 0.00   
];

dx = 0.012;  
dy = 0.0;   

for i = 1:length(x_marks)

    x_target = x_marks(i);
    y_target = y_marks(i);
    x_text = x_target + dx;
    y_text = y_target + dy;
    
    x_target_norm = true_X_norm + (x_target - x_lim(1)) / (x_lim(2) - x_lim(1)) * true_W_norm;
    y_target_norm = true_Y_norm + (y_target - y_lim(1)) / (y_lim(2) - y_lim(1)) * true_H_norm;
    
    x_text_norm = true_X_norm + (x_text - x_lim(1)) / (x_lim(2) - x_lim(1)) * true_W_norm;
    y_text_norm = true_Y_norm + (y_text - y_lim(1)) / (y_lim(2) - y_lim(1)) * true_H_norm;

    annotation('arrow', [x_text_norm, x_target_norm], [y_text_norm, y_target_norm], ...
        'Color', mark_colors(i,:), ...
        'LineWidth', 2, ...
        'HeadStyle', 'vback2', ...
        'HeadWidth', 8, ...
        'HeadLength', 8);

    text(x_text+0.005, y_text, sprintf('(%d)', i), ...
         'FontSize', 16, ...
         'FontWeight', 'bold', ...
         'Color', mark_colors(i,:), ...
         'HorizontalAlignment', 'center', ... 
         'VerticalAlignment', 'middle', ...   
         'BackgroundColor', [1 1 1 0.85], ... 
         'Margin', 1);
end

% exportgraphics(gcf, 'chi12=-4.5.pdf', 'ContentType', 'vector');