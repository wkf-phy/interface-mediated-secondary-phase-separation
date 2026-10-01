format long;

%% 1. Parameters
Nx = 256;
Ny = 256;
L = 1;

chi11 = -7.0; 
chi22 = -7.0; 
chi12 = -4.5;
chi_params = [chi11, chi12, chi22];
kappa_params = [chi11, chi12, chi22];

lr = 0.0045^2;
Det = chi11 * chi22 - chi12^2;
L_sim_tilde = 0.5;
r_min = 1e-6;

c1 = 0.1;
c2_initial = 0.06;

%% 2. Pre-relaxation of the 2D field
fprintf('>Computing the pre-relaxed 2D field...\n');

tol_steady = 1e-4;
phi_steady = get_steady_ternary(c1, c2_initial, Nx, Ny, L, chi_params, kappa_params, tol_steady);

%% 3. Radial profile and mass constraints (PDE -> BVP)
phi1_matrix = phi_steady(:,:,1);
phi2_matrix = phi_steady(:,:,2);

mid_y_idx = round(Ny / 2);
phi1_slice = phi1_matrix(mid_y_idx, :);
phi2_slice = phi2_matrix(mid_y_idx, :);

x_full = linspace(0, L, Nx);
mask = (x_full - L/2) >= 0;

x_plus = x_full(mask) - L/2;
phi1_sub = phi1_slice(mask);
phi2_sub = phi2_slice(mask);

x_plus = x_plus(:)'; 
phi1_sub = phi1_sub(:)'; 
phi2_sub = phi2_sub(:)';

V_total = (L_sim_tilde^2)/2; 
M1_target = trapz(x_plus, phi1_sub .* x_plus); 
M2_target = trapz(x_plus, phi2_sub .* x_plus);

phi1_avg_radial = M1_target / V_total;
phi2_avg_radial = M2_target / V_total;

fprintf('Initial radial means: phi1 = %.5f, phi2 = %.5f\n', ...
    phi1_avg_radial, phi2_avg_radial);

%% 4. Initial guess for the radial BVP
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
curr_p_guess = [mu1_edge; mu2_edge];
curr_x = r_mesh;

solinit = struct('x', curr_x, ...
                 'y', curr_y_guess, ...
                 'parameters', curr_p_guess(:));

%% 5. Continue the radial BVP as the mean client concentration increases
step = 1e-4; 
phi2_avg_list_forward = phi2_avg_radial : step : 0.16;

fprintf('Starting BVP continuation at radial mean phi2 = %.4f\n', ...
    phi2_avg_radial);

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
        
        [warnMsg, ~] = lastwarn();
        if ~isempty(warnMsg)
            error('BVP solver warning: %s', warnMsg);
        end

        phi1_core = res.y(1, 1);
        [~, mid_idx] = min(abs(res.y(1, :) - phi1_core/2.0));
        R_eq = res.x(mid_idx);

        surf_p2 = res.y(3, mid_idx);
        bulk_p2 = res.y(3, end);
        
        actual_phi1_avg = trapz(res.x, res.y(1, :) .* res.x) / V_total;
        actual_phi2_avg = trapz(res.x, res.y(3, :) .* res.x) / V_total;
        
        % Store the converged profile and its measured quantities.
        data_forward(idx_store).phi1_avg = actual_phi1_avg;
        data_forward(idx_store).phi2_avg = actual_phi2_avg;
        data_forward(idx_store).l1 = res.parameters(1);
        data_forward(idx_store).l2 = res.parameters(2);
        data_forward(idx_store).phi2_bulk = bulk_p2;
        data_forward(idx_store).phi1_surf = res.y(1, mid_idx);
        data_forward(idx_store).phi2_surf = surf_p2;
        data_forward(idx_store).R_eq = R_eq;
        data_forward(idx_store).r = res.x;
        data_forward(idx_store).phi1 = res.y(1, :);
        data_forward(idx_store).phi2 = res.y(3, :);
        
        fprintf('Converged: radial mean phi2 = %.4f, R = %.4f\n', ...
            p2_avg, R_eq);

        idx_store = idx_store + 1;
        solinit = struct('x', res.x, 'y', res.y, 'parameters', res.parameters(:));
        
    catch ME
        fprintf('BVP continuation stopped at radial mean phi2 = %.4f\n', ...
            p2_avg);
        fprintf('Reason: %s\n', ME.message);
        break;
    end
end
fprintf('BVP continuation completed: %d solutions obtained.\n', ...
    length(data_forward));

%% 6. Global concentrations and radial BVP instability
fprintf('Converting radial profiles to global concentrations...\n');

n = numel(data_forward);
if n == 0
    error('No converged BVP solutions are available for analysis.');
end
          
A_square = L^2;          
phi2_avg_global_vec = zeros(1, n);
l2_vec = [data_forward.l2];

% Account for the area outside the inscribed circle.
for i = 1:n
    r_d = data_forward(i).r;
    p1 = data_forward(i).phi1;
    p2 = data_forward(i).phi2;

    A_circle = pi * r_d(end)^2;
    A_corner = A_square - A_circle;

    mass2_circle = trapz(r_d, p2 .* r_d) * 2 * pi;

    phi2_bulk = p2(end);
    phi2_avg_global_vec(i) = (mass2_circle + phi2_bulk * A_corner) / A_square;
end

% Find the first positive-to-nonpositive slope crossing.
d_l2_d_phi2 = gradient(l2_vec, phi2_avg_global_vec);
idx = find(d_l2_d_phi2(1:end-1) .* d_l2_d_phi2(2:end) <= 0);

if ~isempty(idx)
    x1 = d_l2_d_phi2(idx); 
    x2 = d_l2_d_phi2(idx+1);
    y1 = phi2_avg_global_vec(idx); 
    y2 = phi2_avg_global_vec(idx+1);

    critical_phi2 = y1 - x1 * (y2 - y1) / (x2 - x1);
    fprintf('Slope-crossing estimate: global mean phi2 = %.5f\n', ...
            critical_phi2);
else
    fprintf('No slope crossing found; last converged global mean phi2 = %.5f\n', ...
            phi2_last_converged_global);
end

%% 7. Plot the radial BVP stability coefficient
fig = figure('Position', [100, 200, 550, 500], 'Color', 'w');
hold on; box on;

plot(phi2_avg_global_vec, d_l2_d_phi2, 'r-', 'LineWidth', 2.0);
yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);

if isfinite(critical_phi2)
    plot(critical_phi2, 0, 'kp', ...
        'MarkerSize', 14, 'MarkerFaceColor', 'y');

    text(critical_phi2 + 0.002, 0.6, ...
        sprintf('$\\bar{\\phi}_{2,c} = %.4f$', critical_phi2), ...
        'Interpreter', 'latex', 'FontSize', 20, ...
        'VerticalAlignment', 'middle', 'FontWeight', 'bold');
end

xlabel('Corrected Global $\bar{\phi}_2$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\partial \mu_2 / \partial \bar{\phi}_2$', 'Interpreter', 'latex', 'FontSize', 20);
title('(b) Susceptibility Analysis', 'Interpreter', 'latex', 'FontSize', 20);

set(gca, 'Layer', 'top', 'FontSize', 20, 'LineWidth', 1.2, 'TickDir', 'in', 'TickLabelInterpreter', 'latex');

xlim([0.06,0.1])
ylim([-5,5])
yticks([-5,0,5])
xticks([0.06,0.08,0.1,0.12])
grid on; set(gca, 'GridColor', [0.95 0.95 0.95], 'LineWidth', 1.0);

% exportgraphics(gcf, 'chi12=-4.5_1D.pdf', 'ContentType', 'vector');

%% Local functions
function [mu1_b, mu2_b] = get_mu_bulk(p1, p2, chi11, chi22, chi12)
%GET_MU_BULK Compute the bulk exchange chemical potentials.

    p1 = max(min(p1, 1-1e-12), 1e-12);
    p2 = max(min(p2, 1-1e-12), 1e-12);
    p3 = max(min(1.0 - p1 - p2, 1-1e-12), 1e-12);
    
    mu1_b = log(p1) - log(p3) + chi11*p1 + chi12*p2;
    mu2_b = log(p2) - log(p3) + chi22*p2 + chi12*p1;
end

function dydr = bvp_ode(r, y, p, chi11, chi22, chi12, lr, Det)
%BVP_ODE Evaluate the radial equilibrium equations and mass integrals.

    p1 = y(1); dp1 = y(2);
    p2 = y(3); dp2 = y(4);
    l1 = p(1); l2 = p(2);
    
    [mu1_b, mu2_b] = get_mu_bulk(p1, p2, chi11, chi22, chi12);
    
    rhs1 = -(mu1_b - l1) / lr;
    rhs2 = -(mu2_b - l2) / lr;
    
    nabla2_p1 = (chi22 * rhs1 - chi12 * rhs2) / Det;
    nabla2_p2 = (-chi12 * rhs1 + chi11 * rhs2) / Det;
    
    dydr = [dp1;
            nabla2_p1 - dp1/r;
            dp2;
            nabla2_p2 - dp2/r;
            p1 * r;
            p2 * r];
end

function res = bvp_bc(ya, yb, p, M1, M2)
%BVP_BC Enforce zero gradients and the prescribed radial masses.

    res = [ya(2);         
           ya(4);         
           ya(5);         
           ya(6);         
           yb(2);         
           yb(4);         
           yb(5) - M1;    
           yb(6) - M2];   
end