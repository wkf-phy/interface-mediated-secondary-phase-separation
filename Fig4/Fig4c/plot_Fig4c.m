load('interfacial_phase_separation_snapshots.mat')

%% Plot morphology time series
% Selected snapshot indices and corresponding time interval
time_indices = 1:length(time_points);

num_rows_imgs = 1; 
num_cols_imgs = length(time_indices); 

%% Component colors
C1 = [244, 100, 209] / 255;
C2 = [22, 195, 80] / 255;  
C3 = [18, 22, 20] / 255;   

%% Figure layout
tile_W = 200; 
tile_H = 200; 

gap_X = 10;   
gap_Y = 0;    

margin_L = 20; 
margin_R = 20; 
margin_T = 60; 
margin_B = 20; 

fig_W = margin_L + num_cols_imgs * tile_W + (num_cols_imgs - 1) * gap_X + margin_R;
fig_H = margin_B + num_rows_imgs * tile_H + (num_rows_imgs - 1) * gap_Y + margin_T;

fig = figure('Color', 'w', 'Position', [200, 200, fig_W, fig_H]); 

%% Plot selected snapshots
for j = 1:num_cols_imgs

    t_idx = time_indices(j);
    actual_time = time_points(t_idx);

    pos_X = (margin_L + (j - 1) * (tile_W + gap_X)) / fig_W;
    pos_Y = margin_B / fig_H; 
    norm_W = tile_W / fig_W;
    norm_H = tile_H / fig_H;

    axes('Position', [pos_X, pos_Y, norm_W, norm_H]);
    
    % Extract the concentration field at the selected time point
    if t_idx <= size(phi_store, 1)
        phi_current = squeeze(phi_store(t_idx, :, :, :));
        
        p1 = phi_current(:,:,1); 
        p2 = phi_current(:,:,2);
        p3 = 1 - p1 - p2; 

        R_ch = p1 * C1(1) + p2 * C2(1) + p3 * C3(1);
        G_ch = p1 * C1(2) + p2 * C2(2) + p3 * C3(2);
        B_ch = p1 * C1(3) + p2 * C2(3) + p3 * C3(3);
        
        phi_rgb = cat(3, R_ch, G_ch, B_ch);
        phi_rgb = max(0, min(1, phi_rgb)); 
        
        imagesc(phi_rgb);
        axis image;

        if actual_time == 0
            title_str = '$t = 0$';
        else
            exponent = floor(log10(actual_time));
            mantissa = actual_time / 10^exponent;
            title_str = sprintf('$t = %g \\times 10^{%d}$', mantissa, exponent);
        end
        title(title_str, 'Interpreter', 'latex', 'FontSize', 25);
    else
        text(0.5, 0.5, 'Index Exceeds Data', 'HorizontalAlignment', 'center', ...
             'Interpreter', 'latex', 'FontSize', 12, 'Color', 'r');
    end

    set(gca, 'XTick', [], 'YTick', []);
    box on
    
end

%% Export
exportgraphics(gcf, 'Morphology_TimeSeries_1x4.png', 'Resolution', 300);
