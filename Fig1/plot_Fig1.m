load("Fig1_dynamics_snapshots.mat") % Simulation data
Nx = 256;
Ny = 256;
L = 1;
x = linspace(0, L, Nx);
y = linspace(0, L, Ny);
dx = L / (Nx - 1);
dy = L / (Ny - 1);
ttotal = 30000;
savedata_time = 1000; 
dtime = 10; 
nstep = ttotal / dtime; 

%% Plot concentration fields and line profiles
disp('Generating concentration maps and line profiles...');

% Component colors
C1 = [244, 100, 209] / 255; % phi1: pink
C2 = [22, 195, 80] / 255; % phi2: green
C3 = [18, 22, 20] / 255; % phi3: dark background

target_times = [0, 4000, 10000, 30000];
num_cols = length(target_times);

figure('Name', 'Time evolution', 'Position', [100, 100, 350*num_cols, 600]);

% Layout in normalized figure coordinates
margin_left = 0.06; 
margin_right = 0.02; 
margin_top = 0.08; 
margin_bottom = 0.2; 
gap_x = -0.25; 
gap_y = 0.04;

ax_width = (1 - margin_left - margin_right - (num_cols - 1) * gap_x) / num_cols;
ax_height = (1 - margin_top - margin_bottom - gap_y) / 2;

% Distance from the selected central grid point
mid_idx_x = round(Nx / 2);
mid_idx_y = round(Ny / 2);
r = x(mid_idx_x:end) - x(mid_idx_x);

for i = 1:num_cols
    t = target_times(i);
    idx = find(time_points == t);
    if idx > size(phi_store, 1) || idx < 1
        warning('No snapshot available for t = %g. Skipping.', t);
        continue;
    end

    phi_current = squeeze(phi_store(idx, :, :, :));
    
    % Top row: concentration maps
    current_left = margin_left + (i - 1) * (ax_width + gap_x);
    bottom_row1 = margin_bottom + ax_height + gap_y;
    
    axes('Position', [current_left, bottom_row1, ax_width, ax_height]);
    
    phi1_img = phi_current(:,:,1);
    phi2_img = phi_current(:,:,2);
    phi3_img = phi_current(:,:,3);
    
    custom_rgb = zeros(Nx, Ny, 3);
    custom_rgb(:,:,1) = phi1_img * C1(1) + phi2_img * C2(1) + phi3_img * C3(1);
    custom_rgb(:,:,2) = phi1_img * C1(2) + phi2_img * C2(2) + phi3_img * C3(2);
    custom_rgb(:,:,3) = phi1_img * C1(3) + phi2_img * C2(3) + phi3_img * C3(3);
    custom_rgb = max(min(custom_rgb, 1), 0);
    
    image(x, y, custom_rgb);
    axis equal tight;
    set(gca, 'YDir', 'normal', 'XTick', [], 'YTick', []);

    title(['$t = ', num2str(t), '$'], 'Interpreter', 'latex', 'FontSize', 23);

    hold on;
    plot([x(mid_idx_x), max(x)], [y(mid_idx_y), y(mid_idx_y)], 'w--', 'LineWidth', 1.5);
    hold off;
    
    
    %  Bottom row: concentration profiles along the selected grid indices
    bottom_row2 = margin_bottom;
    axes('Position', [current_left, bottom_row2, ax_width, ax_height]);
    
    line_phi1 = phi_current(mid_idx_x:end, mid_idx_y, 1);
    line_phi2 = phi_current(mid_idx_x:end, mid_idx_y, 2);

    plot(r, line_phi1, 'Color', C1, 'LineWidth', 1.5); 
    hold on;
    plot(r, line_phi2, 'Color', C2, 'LineWidth', 1.5);
    hold off;
    ylim([0, 1]);
    axis square
    xlim([0, max(r)]);

    set(gca, 'XTick', 0:0.1:0.5, 'YTick', 0:0.2:1, ...
        'FontSize', 16, 'TickLabelInterpreter', 'latex', ...
        'YTickLabel', []);
    grid on;
    
    xlabel('$r$', 'Interpreter', 'latex', 'FontSize', 18);

    % Show the y-axis labels and legend in the first column only
    if i == 1
        ylabel('Volume fraction $\phi$', 'Interpreter', 'latex', 'FontSize', 18);

        lgd = legend('$\phi_1$', '$\phi_2$', 'Location', 'east');
        set(lgd, 'Interpreter', 'latex', 'FontSize', 16);
        set(gca, 'YTickLabel', 0:0.2:1);
    end
end

% Export
% exportgraphics(gcf, 'my_figure.pdf', 'ContentType', 'vector'); 

