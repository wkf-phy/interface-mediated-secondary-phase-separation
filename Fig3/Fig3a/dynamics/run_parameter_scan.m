addpath(genpath('.\ccode'));
format long;

%% Grid and simulation parameters
Nx = 128;   
Ny = 128;
L = 1;
x = linspace(0, L, Nx);
y = linspace(0, L, Ny);
dx = L / (Nx - 1);
dy = L / (Ny - 1);

ttotal = 5000000;
savedata_time = 1000; 
dtime = 10;           
nstep = ttotal / dtime; 

%% Model parameters and Fourier grid
D = 1;
Ldimension = 0.45 * 1e-2;
chi11 = -7;
chi12 = -4.5;   
chi22 = -7;
kappa11 = chi11;
kappa12 = chi12;
kappa22 = chi22;
A = 0.5 * max(max(-[chi11, chi12, chi22]));

% Prepare wave numbers shared by all simulations
[Kx, Ky, k2, k4] = prepare_fft(Nx, Ny, dx, dy);
Kx = Kx(1:end-1);
Ky = Ky(1:end-1);
[kx, ky] = meshgrid(Kx, Ky);

%% Sample parameter sets
c1_range = [0.1, 0.3];
c2_range = [0.06, 0.18];
num_samples = 20; 
rng(1);

[c1_rand, phi2_bar_rand] = generate_lhs_params(c1_range, c2_range, num_samples); 
phase_results = zeros(num_samples, 1);

fprintf('Starting a parallel parameter scan with %d ternary simulations...\n', num_samples);
fprintf('--------------------------------------------------\n');

%% Run simulations in parallel
parfor i_sim = 1:num_samples
    c1 = c1_rand(i_sim);
    phi2_bar = phi2_bar_rand(i_sim);
    c2 = phi2_bar / (1 - c1); 
    
    fprintf('[%d/%d] Starting: c1 = %.3f, phi2_bar = %.3f, c2 = %.3f\n', ...
        i_sim, num_samples, c1, phi2_bar, c2);
        
    % Relax the scaffold droplet before adding the client
    phi1_steady = get_steady_droplet(c1, Nx, Ny, L, 1e-4, 0);
    
    % Replace a fraction of the solvent with the client, including small noise
    phi2_init_field = c2 * ones(Nx, Ny) + phi2_bar / 10 * (rand(Nx, Ny) - 0.5); 
    phi3 = (1 - phi1_steady) - phi2_init_field .* (1 - phi1_steady);
    
    phi1 = phi1_steady;
    phi2 = 1 - phi3 - phi1;
    phi = cat(3, phi1, phi2, phi3);
    
    ttime = 0.0;
    tdraw = savedata_time;
    target_stop_step = inf; 
    phi_old_draw = phi;
    num_phases = 0;
    
    % Time evolution
    for istep = 1:nstep
        ttime = ttime + dtime;
        
        phi1_current = phi(:,:,1);
        phi2_current = phi(:,:,2);
        
        phi1k = fft2(phi1_current);   
        phi2k = fft2(phi2_current);   
        phik = cat(3, phi1k, phi2k);
        
        % Bulk chemical potentials
        mub1 = log(phi1_current) - log(1-phi1_current-phi2_current) + chi11.*phi1_current + chi12.*phi2_current;
        mub1k = fft2(mub1); 
        mub2 = log(phi2_current) - log(1-phi1_current-phi2_current) + chi22.*phi2_current + chi12.*phi1_current;
        mub2k = fft2(mub2); 
            
        % Gradient contributions
        Lphi1k = -k2 .* Ldimension^2 .* (kappa11.*phi1k + kappa12.*phi2k);
        Lphi2k = -k2 .* Ldimension^2 .* (kappa22.*phi2k + kappa12.*phi1k);
        
        muk = cat(3, mub1k + Lphi1k, mub2k + Lphi2k);

        dmudx = real(ifft2(1i .* kx .* muk));  
        dmudy = real(ifft2(1i .* ky .* muk));
        
        dPix = phi(:,:,1).*dmudx(:,:,1) + phi(:,:,2).*dmudx(:,:,2);
        dPiy = phi(:,:,1).*dmudy(:,:,1) + phi(:,:,2).*dmudy(:,:,2);
        
        J1x = phi(:,:,1).*(dmudx(:,:,1) - dPix);
        J1y = phi(:,:,1).*(dmudy(:,:,1) - dPiy);
        J2x = phi(:,:,2).*(dmudx(:,:,2) - dPix);
        J2y = phi(:,:,2).*(dmudy(:,:,2) - dPiy);
        
        dJ1k = 1i.*kx.*fft2(J1x) + 1i.*ky.*fft2(J1y);
        dJ2k = 1i.*kx.*fft2(J2x) + 1i.*ky.*fft2(J2y);
        dJk = cat(3, dJ1k, dJ2k);
        
        % Semi-implicit spectral update.
        Nonlineark = Ldimension^2 .* dJk + A .* Ldimension^4 .* k4 .* phik;
        phik = (Nonlineark.*dtime + phik) ./ (1.0 + A.*k4.*Ldimension^4.*dtime);
        
        phi_temp = real(ifft2(phik));
        phi = cat(3, phi_temp(:,:,1), phi_temp(:,:,2), 1 - phi_temp(:,:,1) - phi_temp(:,:,2));
        
         % Check convergence
        if ttime > (tdraw - dtime/2)
            tdraw = tdraw + savedata_time;
            max_delta_phi = max(abs(phi(:) - phi_old_draw(:)));
            phi_old_draw = phi;
            
            if max_delta_phi < 1e-4
                break;
            end

            if target_stop_step == inf && any(phi(:,:,1) > 0.7, 'all') && any(phi(:,:,2) > 0.7, 'all')
                target_stop_step = istep + 1000; 
            end
        end
        
        if istep >= target_stop_step
            break; 
        end
    end
    
    % Classify the final state by component enrichment.
    has_phase1 = any(phi(:,:,1) > 0.7, 'all');
    has_phase2 = any(phi(:,:,2) > 0.7, 'all');
    has_phase3 = any(phi(:,:,3) > 0.7, 'all');
    num_phases = has_phase1 + has_phase2 + has_phase3;
    
    phase_results(i_sim) = num_phases;
    fprintf('[%d/%d] Finished at step %d; phases detected: %d\n', ...
        i_sim, num_samples, istep, num_phases);
end

%% 5. Plot and save the random-scan results
fprintf('Random scan completed. Plotting the phase classifications...\n');

fig = figure('Name', 'Rapid client addition: parameter scan', ...
    'Position', [100, 100, 600, 500]);
hold on; grid on;

color_1phase = [0, 0.4470, 0.7410];     
color_2phase = [0.4660, 0.6740, 0.1880]; 
color_3phase = [0.9290, 0.6940, 0.1250]; 
sz = 40; 

h = []; labels = {}; 

idx1 = (phase_results == 1);
if any(idx1)
    h(end+1) = scatter(c1_rand(idx1), phi2_bar_rand(idx1), sz, color_1phase, ...
        'filled', 'MarkerEdgeColor', color_1phase, 'LineWidth', 1);
    labels{end+1} = 'One phase';
end

idx2 = (phase_results == 2);
if any(idx2)
    h(end+1) = scatter(c1_rand(idx2), phi2_bar_rand(idx2), sz, color_2phase, ...
        'filled', 'MarkerEdgeColor', color_2phase, 'LineWidth', 1);
    labels{end+1} = 'Two phases';
end

idx3 = (phase_results == 3);
if any(idx3)
    h(end+1) = scatter(c1_rand(idx3), phi2_bar_rand(idx3), sz, color_3phase, ...
        'filled', 'MarkerEdgeColor', color_3phase, 'LineWidth', 1);
    labels{end+1} = 'Three phases';
end

xlabel('Initial mean scaffold concentration, c_1', 'FontSize', 12);
ylabel('Global mean client concentration, \phi_2', 'FontSize', 12);
title('Phase classification after rapid client addition', 'FontSize', 14);

legend(h, labels, 'Location', 'best', 'FontSize', 11);

box on;
hold off;

phi1 = c1_rand;
phi2 = phi2_bar_rand;
phase = phase_results;
save('phase_points.mat', 'phi1', 'phi2', 'phase');