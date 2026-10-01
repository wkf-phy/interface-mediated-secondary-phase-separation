clear; clc;

%% Parameters and component colors
chi11 = -7.0; 
chi22 = -7.0;

num_rows_imgs = 2; 
num_cols_imgs = 2; 

C1 = [244, 100, 209] / 255; % Scaffold
C2 = [22, 195, 80] / 255;   % Client
C3 = [18, 22, 20] / 255;    % Solvent

%% Figure layout

% Tile dimensions, spacing, and margins in pixels
tile_W = 200;   
tile_H = 200;   
gap_X = 5;      
gap_Y = 5;      
margin_L = 90;  
margin_R = 20;  
margin_T = 90;  
margin_B = 20;  

fig_W = margin_L + num_cols_imgs * tile_W + (num_cols_imgs - 1) * gap_X + margin_R;
fig_H = margin_B + num_rows_imgs * tile_H + (num_rows_imgs - 1) * gap_Y + margin_T;

fig = figure('Color', 'w', 'Position', [100, 50, fig_W, fig_H]); 

filename_list = {'case1_no_secondary_phase_separation.mat', ...
    'case2_interface_phase_separation.mat', ...
    'case3_interface_phase_separation.mat', ...
    'case4_bulk_phase_separation.mat'};

% Label colors matching the phase diagram
mark_colors = [
    1.00, 0.00, 0.00;  
    0.00, 0.65, 0.00;  
    0.00, 0.45, 0.85; 
    0.60, 0.20, 0.80   
];

%% Plot morphology maps
for i = 1:num_rows_imgs
    for j = 1:num_cols_imgs

        pos_X = (margin_L + (j - 1) * (tile_W + gap_X)) / fig_W;
        pos_Y = (margin_B + (num_rows_imgs - i) * (tile_H + gap_Y)) / fig_H;
        norm_W = tile_W / fig_W;
        norm_H = tile_H / fig_H;

        axes('Position', [pos_X, pos_Y, norm_W, norm_H]);
        
        filename = filename_list{j+(i-1)*2};
        
        if exist(filename, 'file')
            data = load(filename);
            p1 = data.phi(:,:,1); 
            p2 = data.phi(:,:,2);
            p3 = 1 - p1 - p2; 
            
            % Blend component colors using local volume fractions
            R_ch = p1 * C1(1) + p2 * C2(1) + p3 * C3(1);
            G_ch = p1 * C1(2) + p2 * C2(2) + p3 * C3(2);
            B_ch = p1 * C1(3) + p2 * C2(3) + p3 * C3(3);
            
            phi_rgb = cat(3, R_ch, G_ch, B_ch);
            phi_rgb = max(0, min(1, phi_rgb)); 
            
            imagesc(phi_rgb);
            axis image; 

            hold on;
            idx_mark = j + (i - 1) * 2; 
            xl = xlim; 
            yl = ylim;
            x_pos = xl(1) + 0.02 * (xl(2) - xl(1));
            y_pos = yl(1) + 0.02 * (yl(2) - yl(1));

            text(x_pos, y_pos, sprintf('(%d)', idx_mark), ...
                 'Color', mark_colors(idx_mark, :), ...
                 'FontSize', 14, 'FontWeight', 'bold', ...
                 'BackgroundColor', [1 1 1 0.85], 'Margin', 2, ...
                 'HorizontalAlignment', 'left', ...
                 'VerticalAlignment', 'top');     
        else
            text(0.5, 0.5, 'Data Missing', 'HorizontalAlignment', 'center', ...
                 'Interpreter', 'latex', 'FontSize', 11);
        end
        
        set(gca, 'XTick', [], 'YTick', []);
        box on
    end
end

%% Export
% exportgraphics(gcf, 'Morphology_Diagram.png', 'Resolution', 600);