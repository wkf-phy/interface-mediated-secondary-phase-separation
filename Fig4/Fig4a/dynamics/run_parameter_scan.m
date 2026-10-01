addpath(genpath('.\ccode'));
format long;

%% Grid and simulation parameters
Nx = 128;
Ny = 128;
L = 1;
x = linspace(0, L, Nx);
y = linspace(0, L, Nx);
dx = L / (Nx - 1);
dy = L / (Ny - 1);

ttotal = 5000000;
savedata_time = 1000; 
dtime = 10; 
nstep = ttotal / dtime; 


%% Model parameters and Fourier grid
D_const = 1;
Ldimension = 0.45 * 1e-2;
chi11 = -7;
chi12 = -4.5;
chi22 = -7;
kappa11 = chi11;
kappa12 = chi12;
kappa22 = chi22;
A = 0.5 * max(abs([chi11, chi12, chi22]));

% Prepare wave numbers shared by all simulations
[Kx, Ky, k2, k4] = prepare_fft(Nx, Ny, dx, dy);
Kx = Kx(1:end-1);
Ky = Ky(1:end-1);
[kx, ky] = meshgrid(Kx, Ky);

% Precompute constant coefficients
coeff_L2 = Ldimension^2;
coeff_L4_A = A * Ldimension^4;
denom = 1.0 + coeff_L4_A .* k4 .* dtime; 

%% Sample parameter sets
c1_range = [0.06, 0.25];
c2_range = [0.06, 0.25];
N_samples = 20; 
rng(1);
[c1_list, c2_list] = generate_lhs_params(c1_range, c2_range, N_samples); 

phase_list = zeros(N_samples,1);

D = parallel.pool.DataQueue;
clear printProgress; 
afterEach(D, @(info) printProgress(info, N_samples));

fprintf('Starting a parallel parameter scan with %d ternary simulations...\n', N_samples);
fprintf('--------------------------------------------------\n');

%% Run simulations in parallel
parfor i = 1:N_samples
    t_start = tic; 
    
    c1 = c1_list(i);
    c2 = c2_list(i);
    
     % Initialize concentration fields
    [phi1_init] = initializePhi(Nx, Ny, c1);
    [phi2_init] = initializePhi(Nx, Ny, c2);
    
    phi = zeros(Nx, Ny, 3);
    phi(:,:,1) = phi1_init;
    phi(:,:,2) = phi2_init;
    phi(:,:,3) = 1 - phi(:,:,1) - phi(:,:,2);
    
    phik = complex(zeros(Nx, Ny, 2));
    muk  = complex(zeros(Nx, Ny, 2));
    dmudx = zeros(Nx, Ny, 2);
    dmudy = zeros(Nx, Ny, 2);
    dJk   = complex(zeros(Nx, Ny, 2));
    
    ttime = 0.0;
    tdraw = savedata_time;
    target_stop_step = inf; 
    phi_old_draw = phi;
    
    final_step = nstep; 
    
    % Time evolution
    for istep = 1:nstep
        ttime = ttime + dtime;
        
        phik(:,:,1) = fft2(phi(:,:,1));   
        phik(:,:,2) = fft2(phi(:,:,2));   
        
        % Bulk chemical potentials
        mub1 = log(phi(:,:,1))-log(1-phi(:,:,1)-phi(:,:,2)) + chi11.*phi(:,:,1) + chi12.*phi(:,:,2);
        mub2 = log(phi(:,:,2))-log(1-phi(:,:,1)-phi(:,:,2)) + chi22.*phi(:,:,2) + chi12.*phi(:,:,1);
        
        % Gradient contributions
        Lphi1 = real(ifft2(-k2 .* coeff_L2 .* (kappa11.*phik(:,:,1) + kappa12.*phik(:,:,2))));
        Lphi2 = real(ifft2(-k2 .* coeff_L2 .* (kappa22.*phik(:,:,2) + kappa12.*phik(:,:,1))));
        
        muk(:,:,1) = fft2(mub1) + fft2(Lphi1);
        muk(:,:,2) = fft2(mub2) + fft2(Lphi2);

        dmudx(:,:,1) = real(ifft2(1i*kx.*muk(:,:,1)));  
        dmudx(:,:,2) = real(ifft2(1i*kx.*muk(:,:,2)));
        dmudy(:,:,1) = real(ifft2(1i*ky.*muk(:,:,1)));
        dmudy(:,:,2) = real(ifft2(1i*ky.*muk(:,:,2)));
        
        dPix = phi(:,:,1).*dmudx(:,:,1) + phi(:,:,2).*dmudx(:,:,2);
        dPiy = phi(:,:,1).*dmudy(:,:,1) + phi(:,:,2).*dmudy(:,:,2);
        
        J1x = phi(:,:,1).*(dmudx(:,:,1)-dPix);
        J1y = phi(:,:,1).*(dmudy(:,:,1)-dPiy);
        J2x = phi(:,:,2).*(dmudx(:,:,2)-dPix);
        J2y = phi(:,:,2).*(dmudy(:,:,2)-dPiy);

        dJk(:,:,1) = 1i*kx.*fft2(J1x) + 1i*ky.*fft2(J1y);
        dJk(:,:,2) = 1i*kx.*fft2(J2x) + 1i*ky.*fft2(J2y);
        
        Nonlineark = coeff_L2.*dJk + coeff_L4_A.*k4.*phik;
        
        % Semi-implicit spectral update
        phik = (Nonlineark.*dtime + phik) ./ denom;
        
        phi(:,:,1) = real(ifft2(phik(:,:,1)));
        phi(:,:,2) = real(ifft2(phik(:,:,2)));
        phi(:,:,3) = 1 - phi(:,:,1) - phi(:,:,2);
        
        % Check convergence
        if ttime > (tdraw - dtime/2)
            tdraw = tdraw + savedata_time;
            max_delta_phi = max(abs(phi(:) - phi_old_draw(:)));
            phi_old_draw = phi;
            
            if max_delta_phi < 1e-4
                final_step = istep;
                break;
            end
            
            if target_stop_step == inf && any(phi(:,:,1) > 0.7, 'all') && any(phi(:,:,2) > 0.7, 'all')
                target_stop_step = istep + 1000;
            end
        end
        
        if istep >= target_stop_step
            final_step = istep;
            break; 
        end
    end
    
    % Classify the final state by component enrichment
    has_phase1 = any(phi(:,:,1) > 0.7, 'all');
    has_phase2 = any(phi(:,:,2) > 0.7, 'all');
    has_phase3 = any(phi(:,:,3) > 0.7, 'all');
    
    num_phases = has_phase1 + has_phase2 + has_phase3;
    phase_list(i) = num_phases;
    
    t_elapsed = toc(t_start);
    
    info = struct('c1', c1, 'c2', c2, 'final_step', final_step, ...
                  'num_phases', num_phases, 'time', t_elapsed);
    send(D, info);
end

%% Plot the phase classifications
fprintf('\nParameter scan completed. Plotting the phase classifications...\n');

figure('Name', 'c1-c2 Phase Diagram', 'Position', [100, 100, 600, 500]);
hold on;
grid on;

color_1phase = [0, 0.4470, 0.7410];    
color_2phase = [0.4660, 0.6740, 0.1880]; 
color_3phase = [0.9290, 0.6940, 0.1250]; 

sz = 15; 

h = [];      
labels = {}; 

idx1 = (phase_list == 1);
if any(idx1)
    h(end+1) = scatter(c1_list(idx1), c2_list(idx1), sz, color_1phase, ...
        'filled', 'MarkerEdgeColor', color_1phase, 'LineWidth', 1);
    labels{end+1} = 'One phase';
end

idx2 = (phase_list == 2);
if any(idx2)
    h(end+1) = scatter(c1_list(idx2), c2_list(idx2), sz, color_2phase, ...
        'filled', 'MarkerEdgeColor', color_2phase, 'LineWidth', 1);
    labels{end+1} = 'Two phases';
end

idx3 = (phase_list == 3);
if any(idx3)
    h(end+1) = scatter(c1_list(idx3), c2_list(idx3), sz, color_3phase, ...
        'filled', 'MarkerEdgeColor', color_3phase, 'LineWidth', 1);
    labels{end+1} = 'Three phases';
end

xlabel('Initial concentration, c_1', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Initial concentration, c_2', 'FontSize', 12, 'FontWeight', 'bold');
title('Phase classification from parameter scan', 'FontSize', 14);
legend(h, labels, 'Location', 'best', 'FontSize', 11);

box on;
hold off;

%% Save scan data
fprintf('\nSaving phase-scan data...\n');
phi1 = c1_list;
phi2 = c2_list;
phase = phase_list;
save('phase_points.mat', 'phi1', 'phi2', 'phase');
fprintf('Data saved to: phase_points.mat\n');

%% Local function
function printProgress(info, total_tasks)
    persistent completedTasks;

    if isempty(completedTasks)
        completedTasks = 0;
    end

    completedTasks = completedTasks + 1;
    remainTasks = total_tasks - completedTasks;

     fprintf('[Progress: %3d/%d | Remaining: %3d] Completed [c1=%.2f, c2=%.2f] | Final step: %d | Phases: %d | Time: %.2f s\n', ...
        completedTasks, total_tasks, remainTasks, ...
        info.c1, info.c2, info.final_step, info.num_phases, info.time);
end