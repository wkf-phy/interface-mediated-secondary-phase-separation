function phi1_steady = get_steady_droplet(c1, Nx, Ny, L, tol_steady, plot_flag)
%GET_STEADY_DROPLET Relax a binary mixture to obtain a scaffold droplet.
% c1 is the mean scaffold concentration. Set plot_flag to 0 to disable plots.

    if nargin < 4
        error('At least c1, Nx, Ny, and L must be provided.');
    end
    if nargin < 5 || isempty(tol_steady)
        tol_steady = 1e-4;
    end
    if nargin < 6 || isempty(plot_flag)
        plot_flag = 1;
    end

    Ldimension = 0.45 * 1e-2;
    chi = [0, 3.5; 
           3.5, 0];
    chit = chi(1,1) - chi(1,2) - chi(2,1);
    A = 0.5 * max(max(chi));

    dx = L / (Nx - 1);
    dy = L / (Ny - 1);
    x = linspace(0, L, Nx);
    y = linspace(0, L, Ny);

    dtime = 1;
    max_steps = 50000;
    savedata_time = 1; 
    tdraw = savedata_time;    

    [Kx, Ky, k2, k4] = prepare_fft(Nx, Ny, dx, dy);
    Kx = Kx(1:end-1);
    Ky = Ky(1:end-1);
    [kx, ky] = meshgrid(Kx, Ky);

    phi1 = initializePhi2(Nx, Ny, c1);
    
    if plot_flag
        figure('Name', 'Steady State Generator');
        imagesc(x, y, phi1)
        title('Initial state');
        colorbar;
        axis equal tight;
        drawnow;
    end

    ttime = 0.0;
    
    % Evolve the binary Cahn–Hilliard field.
    for istep = 1:max_steps
        ttime = ttime + dtime;
        phi1_old = phi1; 
        
        phi1k = fft2(phi1);    
        mub1 = log(phi1) - log(1-phi1) + chit .* phi1;
        fhk = fft2(mub1); 
        Lphik = -k2 .* Ldimension^2 .* chit .* phi1k;
        muk = fhk + Lphik;

        dmudx = real(ifft2(1i * kx .* muk));
        dmudy = real(ifft2(1i * ky .* muk));
        
        dPix = phi1 .* dmudx;
        dPiy = phi1 .* dmudy;

        J1x = phi1 .* (dmudx - dPix);
        J1y = phi1 .* (dmudy - dPiy);

        dJ1k = 1i * kx .* fft2(J1x) + 1i * ky .* fft2(J1y);
        Nonlineark = Ldimension^2 .* dJ1k + A * Ldimension^4 .* k4 .* phi1k;
  
        phi1k = (Nonlineark .* dtime + phi1k) ./ (1.0 + A .* k4 .* Ldimension^4 .* dtime);
        phi1 = real(ifft2(phi1k));
        
        % Check the change between consecutive time steps.
        max_diff = max(abs(phi1(:) - phi1_old(:)));
        if max_diff < tol_steady
            fprintf('Binary field converged at t = %.2f (step %d); max change = %.2e\n', ...
                ttime, istep, max_diff);
            if plot_flag
                imagesc(x, y, phi1)
                title(['稳态到达，时间 ', num2str(ttime)]);
                colorbar;
                axis equal tight;
                drawnow;
            end
            break; 
        end

        if plot_flag && (ttime > tdraw)
            tdraw = tdraw + savedata_time;
            imagesc(x, y, phi1)
            title(['Time ', num2str(ttime)]);
            colorbar;
            axis equal tight;
            drawnow; 
        end
    end
    
    if istep == max_steps
        warning('Maximum number of steps (%d) reached before convergence (change: %.2e).', max_steps, max_diff);
    end

    phi1_steady = phi1;
end

function phi = initializePhi2(Nx, Ny, c1)
%INITIALIZEPHI2 Create a droplet profile and match its mean to c1.

    dx = 1 / Nx; 
    x = linspace(0, 1, Nx);
    y = linspace(0, 1, Ny);
    [X, Y] = meshgrid(x, y);
    
    x0 = 0.5;
    y0 = 0.5;
    epsilon = 4 * dx; 
    
    if c1 > 0 && c1 < 1
        radius_est = sqrt(c1 / pi);
    else
        radius_est = 0;
    end
    
    r = sqrt((X - x0).^2 + (Y - y0).^2);
    
    phi_in = 0.95; 
    phi_out = 0.05;
    
    phi_temp = (phi_in + phi_out)/2 + ...
               (phi_in - phi_out)/2 * tanh((radius_est - r) / (sqrt(2) * epsilon));
               
    current_mean = mean(phi_temp(:));
    mass_correction = c1 - current_mean;
    phi = phi_temp + mass_correction;
    
    phi = phi'; 
end