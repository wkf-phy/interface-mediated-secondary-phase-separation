addpath(genpath('.\ccode'));
load('droplet_size_density_analysis.mat')

%% Time points and figure setup
t_list = [1000, 10000, 100000, 150000];

colors = lines(length(t_list)); 
figure('Position', [200, 200, 700, 500], 'Color', 'w');
hold on;
box on;

h_scatter = zeros(length(t_list), 1); 
legend_strs = cell(length(t_list), 1);

%% Analyze droplets at selected time points
for i = 1:length(t_list)

    R_phys = R_list{i};
    rho_phys = rho_list{i};
    
    fprintf('Time t = %6d: %d droplets identified.\n', ...
        t, length(R_phys));

    h_scatter(i) = plot(R_phys, rho_phys, 'o', ...
                    'MarkerSize', 9, ...
                    'MarkerEdgeColor', 'k', ...
                    'MarkerFaceColor', colors(i,:), ...
                    'LineWidth', 1.2, ...
                    'LineStyle', 'none');
                           
    exponent = floor(log10(t));         
    mantissa = t / 10^exponent;  

    if mantissa == 1
        legend_strs{i} = sprintf('$t = 10^{%d}$', exponent);
    else
        legend_strs{i} = sprintf('$t = %g \\times 10^{%d}$', mantissa, exponent);
    end
end

%% Figure formatting
xlabel('$R$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('Linear density $\rho$', 'Interpreter', 'latex', 'FontSize', 20);

set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'latex', 'LineWidth', 1, 'TickDir', 'in');

% Critical linear density
yline(0.0127, '--', 'Color','r', 'LineWidth',2);
text(0.08, 0.0132, {'Critical Linear density $\rho$'},'Interpreter', 'latex','FontSize', 18, ...
    'Color', 'r', 'HorizontalAlignment', 'center','BackgroundColor', 'none','EdgeColor', 'none');  
legend(h_scatter, legend_strs, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best');
hold off;

%% Export
% exportgraphics(gcf, 'Droplets_analyse.pdf', 'ContentType', 'vector');


