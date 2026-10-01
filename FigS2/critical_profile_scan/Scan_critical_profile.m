addpath(genpath('.\ccode'));
format long;
clear; clc; close all;

%% Model and grid parameters
Nx = 256; Ny = 256; L = 1;
chi11 = -7.0; chi22 = -7.0; chi12 = -4.5;
chi_params = [chi11, chi12, chi22];
kappa_params = [chi11, chi12, chi22];
lr = 0.0045^2;
Det = chi11 * chi22 - chi12^2;
L_sim_tilde = 0.5;
r_min = 1e-6;
V_total = 0.5 * (L_sim_tilde^2);
A_square = 1.0^2; 
c2_initial = 0.06; 

%% Scan scaffold concentration
c1_scan_list = 0.07:0.01:0.20; 

critical_points = []; 
critical_profiles = struct([]); 

for idx_c1 = 1:length(c1_scan_list)
    c1_curr = c1_scan_list(idx_c1);
    fprintf('Processing mean scaffold concentration phi1 = %.2f\n', c1_curr);
    
    %% Generate the pre-relaxed 2D concentration field
    phi_steady = get_steady_ternary(c1_curr, c2_initial, Nx, Ny, L, chi_params, kappa_params, 1e-4);
    
    %% Extract radial profiles and mass constraints
    mid_y_idx = round(Ny / 2);
    x_full = linspace(0, L, Nx);
    mask = (x_full - L/2) >= 0;
    
    x_plus = x_full(mask) - L/2; x_plus = x_plus(:)';
    phi1_sub = phi_steady(mid_y_idx, mask, 1); phi1_sub = phi1_sub(:)';
    phi2_sub = phi_steady(mid_y_idx, mask, 2); phi2_sub = phi2_sub(:)';
    
    M1_target = trapz(x_plus, phi1_sub .* x_plus); 
    M2_target = trapz(x_plus, phi2_sub .* x_plus);
    phi1_avg_pde = M1_target / V_total;
    phi2_avg_pde = M2_target / V_total;
    
    %% Construct the initial BVP guess
    r_mesh = linspace(r_min, L_sim_tilde, 2000); 
    dr = r_mesh(2) - r_mesh(1);
    curr_y_guess = zeros(6, length(r_mesh));
    
    curr_y_guess(1, :) = interp1(x_plus, phi1_sub, r_mesh, 'pchip', 'extrap');
    curr_y_guess(3, :) = interp1(x_plus, phi2_sub, r_mesh, 'pchip', 'extrap');
    curr_y_guess(2, :) = gradient(curr_y_guess(1, :), dr);
    curr_y_guess(4, :) = gradient(curr_y_guess(3, :), dr);
    curr_y_guess(5, :) = cumtrapz(r_mesh, curr_y_guess(1, :) .* r_mesh);
    curr_y_guess(6, :) = cumtrapz(r_mesh, curr_y_guess(3, :) .* r_mesh);
    
    [mu1_edge, mu2_edge] = get_mu_bulk(curr_y_guess(1, end), curr_y_guess(3, end), chi11, chi22, chi12);
    solinit = struct('x', r_mesh, 'y', curr_y_guess, 'parameters', [mu1_edge; mu2_edge]);
    
    %% BVP continuation scan
    step = 0.0001; 
    phi2_avg_list_forward = phi2_avg_pde : step : 0.16; 
    
    data_forward = struct([]);
    idx_store = 1;
    options = bvpset('NMax', 100000, 'RelTol', 1e-4);
    
    for p2_avg = phi2_avg_list_forward
        M2_current = p2_avg * V_total;
        fun_curr = @(r, y, p) bvp_ode(r, y, p, chi11, chi22, chi12, lr, Det);
        bc_curr  = @(ya, yb, p) bvp_bc(ya, yb, p, M1_target, M2_current);
        
        try
            lastwarn(''); 
            res = bvp4c(fun_curr, bc_curr, solinit, options);
            
            [warnMsg, warnId] = lastwarn();
            % Treat tolerance warnings as non-converged solutions
            if strcmp(warnId, 'MATLAB:bvp4c:RelTolNotMet') || contains(warnMsg, '满足容差要求') || contains(warnMsg, 'tolerance')
                error('触发非收敛警告');
            end
            
            data_forward(idx_store).phi2_avg = p2_avg;
            data_forward(idx_store).l1 = res.parameters(1); 
            data_forward(idx_store).l2 = res.parameters(2);
            data_forward(idx_store).r = res.x;
            data_forward(idx_store).p1 = res.y(1, :);
            data_forward(idx_store).p2 = res.y(3, :);
            
            idx_store = idx_store + 1;
            solinit = struct('x', res.x, 'y', res.y, 'parameters', res.parameters(:));
            
        catch
            break; 
        end
    end
    
    if isempty(data_forward)
        continue;
    end
    
    %% Convert radial averages to global concentrations
    phi2_glob_vec = zeros(1, length(data_forward));
    phi1_glob_vec = zeros(1, length(data_forward));
    l2_vec = [data_forward.l2];
    
    for i = 1:length(data_forward)
        r_d = data_forward(i).r(:);
        p1 = data_forward(i).p1(:);
        p2 = data_forward(i).p2(:);
        A_corner = A_square - pi*(r_d(end)^2);
        m1_circ = trapz(r_d, p1 .* r_d) * 2 * pi;
        m2_circ = trapz(r_d, p2 .* r_d) * 2 * pi;
        phi1_glob_vec(i) = (m1_circ + p1(end) * A_corner) / A_square;
        phi2_glob_vec(i) = (m2_circ + p2(end) * A_corner) / A_square;
    end
    
    %% Identify the critical point
    d_l2_d_phi2 = gradient(l2_vec, phi2_glob_vec);
    cross_idx = find(d_l2_d_phi2(1:end-1) .* d_l2_d_phi2(2:end) <= 0);
    
    if ~isempty(cross_idx)
        idx_crit = cross_idx(1);
        x1 = d_l2_d_phi2(idx_crit); x2 = d_l2_d_phi2(idx_crit+1);
        y2_crit = phi2_glob_vec(idx_crit) - x1 * (phi2_glob_vec(idx_crit+1) - phi2_glob_vec(idx_crit)) / (x2 - x1);
        y1_crit = phi1_glob_vec(idx_crit) - x1 * (phi1_glob_vec(idx_crit+1) - phi1_glob_vec(idx_crit)) / (x2 - x1);
    else
        idx_crit = length(data_forward);
        y1_crit = phi1_glob_vec(end);
        y2_crit = phi2_glob_vec(end);
    end
    
    %% Store the critical profile
    critical_profiles(idx_c1).c1_base = c1_curr;
    critical_profiles(idx_c1).l1 = data_forward(idx_crit).l1;
    critical_profiles(idx_c1).l2 = data_forward(idx_crit).l2;
    critical_profiles(idx_c1).r = data_forward(idx_crit).r;
    critical_profiles(idx_c1).p1 = data_forward(idx_crit).p1;
    critical_profiles(idx_c1).p2 = data_forward(idx_crit).p2;

    critical_points = [critical_points; y1_crit, y2_crit];
end

%% Plot critical points on the phase diagram
figure('Color', 'w', 'Position', [100, 100, 500, 450]);
plot_cartesian_phase_diagram(chi11, chi12, chi22, [], 0.05);
hold on;
if ~isempty(critical_points)
    plot(critical_points(:, 1), critical_points(:, 2), 'b--', 'LineWidth', 1.5);
    plot(critical_points(:, 1), critical_points(:, 2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k');
end
title('Critical Points in the Phase Diagram', ...
    'FontSize', 14);

%% Plot critical concentration profiles
figure('Color', 'w', 'Position', [650, 100, 900, 400]);
color_map = lines(length(c1_scan_list));

% Scaffold profiles
subplot(1, 2, 1); hold on; box on;
for i = 1:length(c1_scan_list)
    if length(critical_profiles) >= i && ~isempty(critical_profiles(i).r)
        plot(critical_profiles(i).r, critical_profiles(i).p1, ...
            'LineWidth', 2, 'Color', color_map(i,:), ...
            'DisplayName', sprintf('\\phi_{1,avg} = %.2f', critical_profiles(i).c1_base));
    end
end
xlabel('Radial distance $r$', ...
    'Interpreter', 'latex', ...
    'FontSize', 12);

ylabel('Scaffold concentration $\phi_1(r)$', ...
    'Interpreter', 'latex', ...
    'FontSize', 12);

title('Critical $\phi_1$ Profiles', ...
    'Interpreter', 'latex', ...
    'FontSize', 14);
legend('Location', 'best', 'FontSize', 11);

% Client profiles
subplot(1, 2, 2); hold on; box on;
for i = 1:length(c1_scan_list)
    if length(critical_profiles) >= i && ~isempty(critical_profiles(i).r)
        plot(critical_profiles(i).r, critical_profiles(i).p2, ...
            'LineWidth', 2, 'Color', color_map(i,:), ...
            'DisplayName', sprintf('\\phi_{1,avg} = %.2f', critical_profiles(i).c1_base));
    end
end
xlabel('Radial distance $r$', ...
    'Interpreter', 'latex', ...
    'FontSize', 12);

ylabel('Client concentration $\phi_2(r)$', ...
    'Interpreter', 'latex', ...
    'FontSize', 12);

title('Critical $\phi_2$ Profiles', ...
    'Interpreter', 'latex', ...
    'FontSize', 14);
legend('Location', 'best', 'FontSize', 11);

fprintf('Critical profiles extracted and plotted successfully.\n');

%% Save critical profiles
save('critical_profile.mat', 'critical_profiles', 'critical_points');
fprintf('Critical profiles saved to: critical_profile.mat\n');

%% Local functions
function [mu1_b, mu2_b] = get_mu_bulk(p1, p2, chi11, chi22, chi12)
    p1 = max(min(p1, 1-1e-12), 1e-12);
    p2 = max(min(p2, 1-1e-12), 1e-12);
    p3 = max(min(1.0 - p1 - p2, 1-1e-12), 1e-12);
    mu1_b = log(p1) - log(p3) + chi11*p1 + chi12*p2;
    mu2_b = log(p2) - log(p3) + chi22*p2 + chi12*p1;
end
function dydr = bvp_ode(r, y, p, chi11, chi22, chi12, lr, Det)
    p1 = y(1); dp1 = y(2); p2 = y(3); dp2 = y(4);
    l1 = p(1); l2 = p(2);
    [mu1_b, mu2_b] = get_mu_bulk(p1, p2, chi11, chi22, chi12);
    rhs1 = -(mu1_b - l1) / lr; rhs2 = -(mu2_b - l2) / lr;
    nabla2_p1 = (chi22 * rhs1 - chi12 * rhs2) / Det;
    nabla2_p2 = (-chi12 * rhs1 + chi11 * rhs2) / Det;
    dydr = [dp1; nabla2_p1 - dp1/r; dp2; nabla2_p2 - dp2/r; p1*r; p2*r];
end
function res = bvp_bc(ya, yb, p, M1, M2)
    res = [ya(2); ya(4); ya(5); ya(6); yb(2); yb(4); yb(5) - M1; yb(6) - M2];   
end