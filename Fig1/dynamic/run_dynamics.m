addpath(genpath('.\ccode'));

format long;

%% Grid and simulation parameters
Nx = 256;
Ny = 256;
L = 1;
x = linspace(0, L, Nx);
y = linspace(0, L, Ny);
dx = L / (Nx - 1);
dy = L / (Ny - 1);

ttotal = 30000;
savedata_time = 1000;    % Interval between saved snapshots
dtime = 10;              % Time step
nstep = ttotal / dtime;  

steps_per_save = round(savedata_time / dtime);       
total_save_points = floor(nstep / steps_per_save) + 1; % Include t = 0

ttime = 0.0;           
flag = 2;

% Model parameters
D = 1;
Ldimension = 0.45 * 1e-2;

chi11 = -7;
chi12 = -4.5;
chi22 = -7;

kappa11 = chi11;
kappa12 = chi12;
kappa22 = chi22;

A = 0.5 * max(abs([chi11, chi12, chi22]));

%% Fourier grid and initial concentration fields
[Kx, Ky, k2, k4] = prepare_fft(Nx, Ny, dx, dy);
Kx = Kx(1:end-1);
Ky = Ky(1:end-1);
[kx, ky] = meshgrid(Kx, Ky);

% Load the binary equilibrium state
base_data = load("phi1_phi3_equilibrium.mat");
phi_base = base_data.phi;
phi = imresize(phi_base, [Nx, Ny], 'bicubic');

% Add the client by replacing a fraction of the solvent
phi2_bar = 0.08;
c2 = phi2_bar / mean(mean(phi(:,:,2)));
phi2 = initializePhi(Nx, Ny, c2);

phi3 = phi(:,:,2) - phi2 .* (phi(:,:,2));
phi1 = phi(:,:,1);
phi2 = 1 - phi3 - phi1;

phi = cat(3, phi1, phi2, phi3); 

%% Allocate snapshot storage
% Array dimensions: snapshot, grid index 1, grid index 2, component
phi_store = zeros([total_save_points, Nx, Ny, 3]);
save_count = 1; 
phi_store(save_count, :, :, :) = phi; 

disp('Initial state stored at t = 0.');

phi2_old = phi(:,:,2);

%% Time evolution
for istep = 1:nstep
    ttime = ttime + dtime;
    
    phi1 = phi(:,:,1);
    phi2 = phi(:,:,2);

    phi1k = fft2(phi1);   
    phi2k = fft2(phi2);   
    phik = cat(3, phi1k, phi2k);

    % Bulk contributions to the exchange chemical potentials
    mub1 = log(phi1) - log(1-phi1-phi2) + chi11.*phi1 + chi12.*phi2;
    mub1k = fft2(mub1); 
    mub2 = log(phi2) - log(1-phi1-phi2) + chi22.*phi2 + chi12.*phi1;
    mub2k = fft2(mub2); 
        
    % Gradient contributions to the exchange chemical potentials
    Lphi1k = -k2 .* Ldimension^2 .* (kappa11.*phi1k + kappa12.*phi2k);
    Lphi1 = real(ifft2(Lphi1k));
    Lphi2k = -k2 .* Ldimension^2 .* (kappa22.*phi2k + kappa12.*phi1k);
    Lphi2 = real(ifft2(Lphi2k));

    mu1k = mub1k + Lphi1k;
    mu2k = mub2k + Lphi2k;
    muk = cat(3, mu1k, mu2k);

    % Chemical-potential gradients
    dmudx = real(ifft2(1i * kx .* muk));  
    dmudy = real(ifft2(1i * ky .* muk));
    
    dPix = phi(:,:,1).*dmudx(:,:,1) + phi(:,:,2).*dmudx(:,:,2);
    dPiy = phi(:,:,1).*dmudy(:,:,1) + phi(:,:,2).*dmudy(:,:,2);

    J1x = phi(:,:,1).*(dmudx(:,:,1) - dPix);
    J1y = phi(:,:,1).*(dmudy(:,:,1) - dPiy);

    J2x = phi(:,:,2).*(dmudx(:,:,2) - dPix);
    J2y = phi(:,:,2).*(dmudy(:,:,2) - dPiy);

    % Divergences in Fourier space
    dJ1k = 1i*kx.*fft2(J1x) + 1i*ky.*fft2(J1y);
    dJ2k = 1i*kx.*fft2(J2x) + 1i*ky.*fft2(J2y);
    dJk = cat(3, dJ1k, dJ2k);

    Nonlineark = Ldimension^2.*dJk + A*Ldimension^4.*k4.*phik;

    % Semi-implicit update in Fourier space
    phik = (Nonlineark.*dtime + phik) ./ (1.0 + A.*k4.*Ldimension^4.*dtime);

    phi_temp = real(ifft2(phik));
    phi_temp = max(min(phi_temp, 1), 0); 
    phi = cat(3, phi_temp(:,:,1), phi_temp(:,:,2), 1 - phi_temp(:,:,1) - phi_temp(:,:,2));

    % Store snapshots at the specified interval
    if mod(istep, steps_per_save) == 0
        save_count = save_count + 1;           
        phi_store(save_count, :, :, :) = phi;  
        fprintf('Stored snapshot %d at t = %g.\n', save_count, ttime);
    end
end

save('Fig1_data.mat', 'phi_store', '-v7.3');
