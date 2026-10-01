function phi_steady = get_steady_ternary(c1, c2, Nx, Ny, L, chi_params, kappa_params, tol_steady)
%GET_STEADY_TERNARY Relax a ternary mixture to a steady concentration field.

    if nargin < 8 || isempty(tol_steady)
        tol_steady = 1e-4;
    end
    
    max_steps = 100000; 
    dtime = 50;         
    check_time = 1000;  

    chi11 = chi_params(1); chi12 = chi_params(2); chi22 = chi_params(3);
    kappa11 = kappa_params(1); kappa12 = kappa_params(2); kappa22 = kappa_params(3);

    dx = L / (Nx - 1);
    dy = L / (Ny - 1);
    Ldimension = 0.45 * 1e-2;
    A = 0.5 * max(abs([chi11, chi12, chi22]));

    [Kx, Ky, k2, k4] = prepare_fft(Nx, Ny, dx, dy);
    Kx = Kx(1:end-1); Ky = Ky(1:end-1);
    [kx, ky] = meshgrid(Kx, Ky);

    % Initialize from a relaxed scaffold droplet.
    fprintf('Preparing the ternary steady state...\n');
    fprintf('Step 1/2: Relaxing the scaffold droplet...\n');
    phi1_steady = get_steady_droplet(c1, Nx, Ny, L, 1e-4, 0);

    phi2_bar = c2; 
    c2_scaled = phi2_bar / mean(mean(1 - phi1_steady));
    
    phi2_init_field = c2_scaled * ones(Nx, Ny) + 1e-3 * (rand(Nx, Ny) - 0.5);
    phi3 = (1 - phi1_steady) - phi2_init_field .* (1 - phi1_steady);
    phi1 = phi1_steady;
    phi2 = 1 - phi3 - phi1;
    phi = cat(3, phi1, phi2, phi3);

    fprintf('Step 2/2: Relaxing the ternary field (check interval: %d; tolerance: %.2e)...\n', ...
        check_time, tol_steady);

    ttime = 0;
    next_check_time = check_time;
    phi_old = phi; % 
    
    for istep = 1:max_steps
        ttime = ttime + dtime;

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

         % Check the maximum field change over one check interval.
        if ttime >= next_check_time
            max_diff = max(abs(phi(:) - phi_old(:)));
            if max_diff < tol_steady
                fprintf('Ternary field converged at t = %.2f (step %d); max change = %.2e\n', ...
                    ttime, istep, max_diff);
                break;
            end
            phi_old = phi; 
            next_check_time = next_check_time + check_time;
        end
    end

    if istep == max_steps
        warning('Maximum number of steps (%d) reached before convergence (change: %.2e).', ...
            max_steps, max_diff);
    end

    phi_steady = phi;
end