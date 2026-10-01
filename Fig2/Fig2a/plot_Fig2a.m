% Plot phase boundaries and simulation results
addpath(genpath('.\selected_simulation_data'));
load('instability_boundary.mat')

% Extract valid boundary data
valid_mask = ~isnan(critical_phi2_list);
x_v = chi12_list(valid_mask);
y_v = critical_phi2_list(valid_mask);

y_min = 0.06; 
y_max = 0.2;
x_min = -4.8;
x_max = 0; 

% Interpolate the boundary for background shading
x_dense = linspace(x_min, x_max, 1000);
y_dense = interp1(x_v, y_v, x_dense, 'linear', 'extrap');

% Background colors
color_reg1 = [0.85, 0.90, 1.00]; 
color_reg2 = [1.00, 0.98, 0.80]; 
color_reg3 = [1.00, 0.90, 0.85]; 
color_reg4 = [0.85, 1.00, 0.85]; 

figure('Position', [150, 150, 800, 600], 'Color', 'w');
hold on; 

% Region 1: x < -4.8
if x_min < -4.8
    fill([x_min, -4.8, -4.8, x_min], [y_min, y_min, y_max, y_max], ...
         color_reg1, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
end

% Region 2: below the boundary
m2 = x_dense >= -4.8 & x_dense <= 0;
x2 = x_dense(m2); y2 = y_dense(m2);
if ~isempty(x2)
    fill([x2, fliplr(x2)], [y2, repmat(y_min, 1, length(y2))], ...
         color_reg2, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
end

% Region 3: above the boundary, -4.8 <= x <= -1.3
m3 = x_dense >= -4.8 & x_dense <= -1.3;
x3 = x_dense(m3); y3 = y_dense(m3);
if ~isempty(x3)
    fill([x3, fliplr(x3)], [y3, repmat(y_max, 1, length(y3))], ...
         color_reg3, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
end

% Region 4: above the boundary, -1.3 <= x <= 0
m4 = x_dense >= -1.3 & x_dense <= 0;
x4 = x_dense(m4); y4 = y_dense(m4);
if ~isempty(x4)
    fill([x4, fliplr(x4)], [y4, repmat(y_max, 1, length(y4))], ...
         color_reg3, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
end

% Overlay boundary points and reference lines
plot(x_v(x_v>-4.9), y_v(x_v>-4.9), '-o', 'LineWidth', 2, 'MarkerSize', 4, ...
    'MarkerFaceColor', 'b', 'Color', 'b'); %不画两相区
yline(0.173, '--', 'LineWidth', 2,'Color','r')
xline(-1.3, '--', 'LineWidth', 2,'Color','b')

box on; grid on;
set(gca, 'Layer', 'top'); 
axis([x_min x_max y_min y_max]); 
set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex', 'GridColor', [0.8 0.8 0.8]);

% Simulation grid and classification thresholds
chi12_array = -5:0.25:0;          
phi2_bar_array = 0.06:0.01:0.2; 
chi11 = -7.0; 
chi22 = -7.0;
PS_threshold = 0.5;       
noise_area_limit = 5;      

color_ps = [0.85, 0.15, 0.15];     % Phase separation detected
color_nops = [1.00, 0.60, 0.00];   % No phase separation detected

num_chi12 = length(chi12_array);
num_phi2_bar = length(phi2_bar_array);

y_crit_interp_scatter = interp1(x_v, y_v, chi12_array, 'linear', 'extrap');
plot_mask = false(num_phi2_bar, num_chi12);

for j = 1:num_chi12
    yc = y_crit_interp_scatter(j); 
    idx_above = find(phi2_bar_array >= yc, 1, 'first');
    
    if isempty(idx_above)
        target_indices = [num_phi2_bar-1, num_phi2_bar];
    elseif idx_above == 1
        target_indices = [1, 2];
    else
        k = idx_above - 1; 
        target_indices = [k-1, k, k+1, k+2]; 
    end
    
    target_indices = target_indices(target_indices >= 1 & target_indices <= num_phi2_bar);
    plot_mask(target_indices, j) = true;
end

ax1 = gca; 

%% Plot simulation results

load('simulation_classification.mat');

for i = 1:length(phase_points)

    chi12 = chi12_points(i);
    phi2_bar = phi2_points(i);

    if phase_points(i)
        scatter(ax1, chi12, phi2_bar, 200, '^', ...
            'MarkerFaceColor', color_ps, ...
            'MarkerEdgeColor', [0.2 0.2 0.2], ...
            'LineWidth', 0.8);
    else
        scatter(ax1, chi12, phi2_bar, 200, 's', ...
            'MarkerFaceColor', color_nops, ...
            'MarkerEdgeColor', [0.2 0.2 0.2], ...
            'LineWidth', 0.8);
    end
end

% Create legend handles
h_ps = plot(ax1, NaN, NaN, '^', 'MarkerSize', 13, 'MarkerFaceColor', color_ps, ...
    'MarkerEdgeColor', [0.2 0.2 0.2], 'LineStyle', 'none');
h_nops = plot(ax1, NaN, NaN, 's', 'MarkerSize', 14, 'MarkerFaceColor', color_nops, ...
    'MarkerEdgeColor', [0.2 0.2 0.2], 'LineStyle', 'none');

text(-3.5, 0.14, 'Three-phase','Interpreter', 'latex','FontSize', 18, ...
    'Color', 'k', 'HorizontalAlignment', 'center','BackgroundColor', 'none','EdgeColor', 'none');  
text(-0.7, 0.12, 'Two-phase','Interpreter', 'latex','FontSize', 18, ...
    'Color', 'k', 'HorizontalAlignment', 'center','BackgroundColor', 'none','EdgeColor', 'none');  
text(-3, 0.183, {'Instability threshold without $\phi_1$'},'Interpreter', 'latex','FontSize', 18, ...
    'Color', 'r', 'HorizontalAlignment', 'center','BackgroundColor', 'none','EdgeColor', 'none');  

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 25);

% Finalize the layout
box on;
lgd = legend(ax1, [h_ps, h_nops], {'Three-phase', 'Two-phase'}, ...
    'Location', 'southeast', 'NumColumns', 1, 'Interpreter', 'latex', 'FontSize', 25);
xlabel('Interaction parameter $\chi_{12}$', 'Interpreter', 'latex', 'FontSize', 30);
ylabel('Volume fraction $\bar{\phi_2}$', 'Interpreter', 'latex', 'FontSize', 30);

ht=title(['Phase diagram at $\bar{\phi}_1 = ', num2str(c1), '$'], ...
          'Interpreter', 'latex', 'FontSize', 30);
set(ht, 'Units', 'normalized', 'Position', [0.5, 1.03, 0]);

pos = get(gca, 'Position'); 
pos(4) = pos(4) * 0.9; 
set(gca, 'Position', pos);
xticks(-5:1:0); 
yticks(0.05:0.05:0.2); 
drawnow; 

x_marks = [-1, -4];
y_marks = [0.17, 0.07];

mark_colors = [
    0.60, 0.20, 0.80  
    0.00, 0.80, 0.00 
];

dx = 0.0;  
dy = 0.007;   

ax_pos = ax1.Position; 
x_lim = ax1.XLim;
y_lim = ax1.YLim;

for i = 1:length(x_marks)
    x_target = x_marks(i);
    y_target = y_marks(i);
    
    x_text = x_target + dx;
    y_text = y_target + dy;
    
    x_target_norm = ax_pos(1) + (x_target - x_lim(1)) / (x_lim(2) - x_lim(1)) * ax_pos(3);
    y_target_norm = ax_pos(2) + (y_target - y_lim(1)) / (y_lim(2) - y_lim(1)) * ax_pos(4);
    
    x_text_norm = ax_pos(1) + (x_text - x_lim(1)) / (x_lim(2) - x_lim(1)) * ax_pos(3);
    y_text_norm = ax_pos(2) + (y_text - y_lim(1)) / (y_lim(2) - y_lim(1)) * ax_pos(4);

    annotation('arrow', [x_text_norm, x_target_norm], [y_text_norm, y_target_norm], ...
        'Color', mark_colors(i,:), ...
        'LineWidth', 2, ...
        'HeadStyle', 'vback2', ...
        'HeadWidth', 8, ...
        'HeadLength', 8);
    
    text(x_text, y_text+0.003, sprintf('(%d)', i), ...
         'FontSize', 16, ...
         'FontWeight', 'bold', ...
         'Color', mark_colors(i,:), ...
         'HorizontalAlignment', 'center', ... 
         'VerticalAlignment', 'middle', ...   
         'BackgroundColor', [1 1 1 0.85], ... 
         'Margin', 1);
end

% exportgraphics(gcf, 'my_figure.pdf', 'ContentType', 'vector');

