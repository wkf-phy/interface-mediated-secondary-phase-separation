addpath(genpath('.\ccode'));
format long;

%% Grid and time parameters
Nx = 256;
Ny = 256;
L = 1;
x = linspace(0, L, Nx);
y = linspace(0, L, Ny);
dx = L / (Nx - 1);
dy = L / (Ny - 1);
ttotal = 5000000;
savedata_time = 1000; 
dtime = 10;          
tdraw = savedata_time; 
ttime = 0.0;        
nstep = ttotal / dtime; 
target_stop_step = inf; 

%% Model parameters
D = 1;
Ldimension = 0.45 * 1e-2;
chi11 = -7;
chi12 = -4.5;
chi22 = -7;
kappa11 = chi11;
kappa12 = chi12;
kappa22 = chi22;

% No secondary phase separation
% c1 = 0.1; 
% phi2_bar = 0.08;

% Interfacial phase separation
c1 = 0.147694; 
phi2_bar = 0.0951886;

% Interfacial phase separation2
% c1 = 0.16282; 
% phi2_bar = 0.120397;

% Bulk phase separation
% c1 = 0.120369; 
% phi2_bar = 0.20321;

c2 = phi2_bar/(1-c1);
A = 0.5 * max(max(-[chi11, chi12, chi22]));

[Kx, Ky, k2, k4] = prepare_fft(Nx, Ny, dx, dy);
Kx = Kx(1:end-1);
Ky = Ky(1:end-1);
[kx, ky] = meshgrid(Kx, Ky);

%% Initialize concentration fields
fprintf('Generating the pre-relaxed scaffold droplet...\n');
phi1_steady = get_steady_droplet(c1, Nx, Ny, L, 1e-4, 0);
phi2_init_field = c2 * ones(Nx, Ny) + 1e-3 * (rand(Nx, Ny) - 0.5); 
phi3 = (1 - phi1_steady) - phi2_init_field .* (1 - phi1_steady);
phi1 = phi1_steady;
phi2 = 1 - phi3 - phi1;
phi = cat(3, phi1, phi2, phi3); 
phi_store = zeros([ttotal/savedata_time, Nx, Ny, 3]);
phi_old_draw = phi;

%% Initialize the figure
figure(1);
set(gcf, 'Position', [100, 100, 1000, 800], 'Color', 'w'); 
tlo = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

%% Cahn-Hilliard dynamics
for istep = 1:nstep
    ttime = ttime + dtime;  
    
    phi2_old = phi(:,:,2);
    
    phi1 = phi(:,:,1);
    phi2 = phi(:,:,2);
    
    phi1k = fft2(phi1);   
    phi2k = fft2(phi2);   
    phik = cat(3, phi1k, phi2k);
    
    mub1 = log(phi1) - log(1-phi1-phi2) + chi11.*phi1 + chi12.*phi2;
    mub1k = fft2(mub1); 
    mub2 = log(phi2) - log(1-phi1-phi2) + chi22.*phi2 + chi12.*phi1;
    mub2k = fft2(mub2); 
        
    Lphi1k = -k2 .* Ldimension^2 .* (kappa11.*phi1k + kappa12.*phi2k);
    Lphi2k = -k2 .* Ldimension^2 .* (kappa22.*phi2k + kappa12.*phi1k);
    
    mu1k = mub1k + Lphi1k;
    mu2k = mub2k + Lphi2k;
    muk = cat(3, mu1k, mu2k);
    
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
    
    Nonlineark = Ldimension^2 .* dJk + A .* Ldimension^4 .* k4 .* phik;
    phik = (Nonlineark.*dtime + phik) ./ (1.0 + A.*k4.*Ldimension^4.*dtime);
    
    phi_temp = real(ifft2(phik));
    phi = cat(3, phi_temp(:,:,1), phi_temp(:,:,2), 1 - phi_temp(:,:,1) - phi_temp(:,:,2));
    
    % Maximum client concentration change over one time step
    max_diff_phi2 = max(abs(phi(:,:,2) - phi2_old), [], 'all');
    
    
    % Store snapshots and update the figure
    if ttime > (tdraw - dtime/2)
        tdraw = tdraw + savedata_time;
        num = round(ttime / savedata_time);
        phi_store(num,:,:,:) = phi;
        
        % Scaffold concentration map
        nexttile(tlo, 1);
        imagesc(x, y, phi(:,:,1)); 
        colorbar; shading interp; axis equal tight;
        title(['Scaffold volume fraction \phi_1, t=', num2str(ttime)]);
        
        % Client concentration map
        nexttile(tlo, 2);
        imagesc(x, y, phi(:,:,2));
        colorbar; shading interp; axis equal tight;
        title(['Client volume fraction \phi_2, t=', num2str(ttime)]);
        subtitle(sprintf('Maximum single-step change: %.2e', max_diff_phi2), ...
            'Color', 'r', 'FontWeight', 'bold');
        
        % Concentration profiles along the central grid row
        nexttile(tlo, 3, [1, 2]);
        mid_y = round(Ny / 2); 
        
        plot(x, phi(mid_y, :, 1), 'r-', 'LineWidth', 1.5); hold on;
        plot(x, phi(mid_y, :, 2), 'b-', 'LineWidth', 1.5);
        plot(x, phi(mid_y, :, 3), 'k--', 'LineWidth', 1.5); hold off;

        grid on;
        axis([0 L 0 1]); 
        xlabel('x'); ylabel('Concentration');
        legend('\phi_1', '\phi_2', '\phi_3', 'Location', 'northeast');
        title(['Central line profiles, t=', num2str(ttime)]);
        
        pause(0.01);
        drawnow; 

        max_delta_phi = max(abs(phi(:) - phi_old_draw(:)));
        phi_old_draw = phi;

        if max_delta_phi < 1e-4
            break;
        end

        if target_stop_step == inf && any(phi(:,:,1) > 0.75, 'all') && any(phi(:,:,2) > 0.75, 'all')
            target_stop_step = istep*2; 
            i
        end
    end

    if istep >= target_stop_step
        break; 
    end
end

%% Save the final state
case_id = 2; % Set to 1, 2, 3, or 4 for the current simulation
case_files = {
    'case1_no_secondary_phase_separation.mat'
    'case2_interface_phase_separation.mat'
    'case3_interface_phase_separation.mat'
    'case4_bulk_phase_separation.mat'
};
save(case_files{case_id}, 'phi', 'ttime');