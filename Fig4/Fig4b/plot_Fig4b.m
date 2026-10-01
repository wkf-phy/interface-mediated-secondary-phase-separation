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

filename_list={'case1_no_phase_separation.mat',...
    'case2_interfacial_adsorption.mat',...
    'case3_interfacial_phase_separation.mat',...
    'case4_interfacial_phase_separation.mat'};

mark_colors = [
    1.00, 0.00, 0.00; 
    0.00, 0.80, 0.00;  
    0.00, 0.00, 1.00;  
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

        filename = filename_list{2*(i-1)+j};
        
        if exist(filename, 'file')
            data = load(filename);
            p1 = data.phi(:,:,1); 
            p2 = data.phi(:,:,2);
            p3 = 1 - p1 - p2; 
            
            R_ch = p1 * C1(1) + p2 * C2(1) + p3 * C3(1);
            G_ch = p1 * C1(2) + p2 * C2(2) + p3 * C3(2);
            B_ch = p1 * C1(3) + p2 * C2(3) + p3 * C3(3);
            
            phi_rgb = cat(3, R_ch, G_ch, B_ch);
            phi_rgb = max(0, min(1, phi_rgb)); 
            
            imagesc(phi_rgb);
            axis image; 

            hold on;

            idx_mark = 2*(i-1)+j; 
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

        if idx_mark == 2
            % Add a one-dimensional concentration profile
            mid_row = round(size(p2, 1) / 2) + 5; 
            plot([70, 90], [mid_row, mid_row], 'w--', 'LineWidth', 1.5);
            
            inset_W = norm_W * 0.4; 
            inset_H = norm_H * 0.3; 
            inset_X = pos_X + norm_W - inset_W - 0.02; 
            inset_Y = pos_Y + norm_H - inset_H - 0.02; 
            ax_inset = axes('Position', [inset_X, inset_Y, inset_W, inset_H]);

            plot(p1(mid_row, 70:94), 'Color', C1, 'LineWidth', 1.5); hold on;
            plot(p2(mid_row, 70:94), 'Color', C2, 'LineWidth', 1.5);

            lgd = legend('$\phi_1$', '$\phi_2$', 'Interpreter', 'latex', 'Location', 'best');
            lgd.FontSize = 8; 
            lgd.ItemTokenSize = [12, 18]; 
            lgd.TextColor = 'w'; 
            legend boxoff;      

            xlim([1, 94-70]);
            ylim([0, 0.85]);
            set(ax_inset, 'Color', [0.25, 0.25, 0.25], ... 
                          'XColor', 'w', ...               
                          'YColor', 'w', ...               
                          'XTick', [], 'YTick', []);      
            box on; 

            uistack(ax_inset, 'top');
            ax_main = findall(gcf, 'type', 'axes', 'Position', [pos_X, pos_Y, norm_W, norm_H]);
            set(gcf, 'CurrentAxes', ax_main(1)); 
        end
    end
end

%% Export
% exportgraphics(gcf, 'Morphology_Diagram.pdf', 'ContentType', 'vector');