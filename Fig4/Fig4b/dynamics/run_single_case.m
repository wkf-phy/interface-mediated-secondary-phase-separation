addpath(genpath('.\ccode'));
format long;

%% Grid and time parameters
Nx = 128;
Ny = 128;
L=1;
x = linspace(0,L,Nx);
y = linspace(0,L,Nx);
dx=L/(Nx-1);
dy=L/(Ny-1);

ttotal=250000;
savedata_time=1000; 
dtime = 10;

tdraw = savedata_time;   
ttime = 0.0;   
nstep = ttotal/dtime;   

%% Model parameters
D=1;
Ldimension=0.45*1e-2;
chi11= -7;
chi12= -4.5;
chi22= -7;

kappa11=chi11;
kappa12=chi12;
kappa22=chi22;

A=0.5*max(abs([chi11,chi12,chi22]));

% Prepare wave numbers
[Kx, Ky, k2, k4] = prepare_fft(Nx, Ny, dx, dy);
Kx = Kx(1:end-1);
Ky = Ky(1:end-1);
[kx, ky] = meshgrid(Kx, Ky);

% No phase separation
% c1 = 0.105873;
% c2 = 0.095714;

% Interfacial adsorption
% c1 = 0.160342;
% c2 = 0.077886;

% Interfacial phase separation
% c1 = 0.179262;
% c2 = 0.111275;

% Interfacial phase separation 2
c1 = 0.107518;
c2 = 0.165370;

%% Initialize concentration fields
[phi1] = initializePhi(Nx, Ny, c1);
[phi2] = initializePhi(Nx, Ny, c2);
phi3=1-phi1-phi2;

phi = cat(3, phi1, phi2,phi3);
phi_store=zeros([ttotal/savedata_time,Nx,Ny,3]);
phi2_old=0;

%% Initialize the figure
figure('Color', 'w', 'Position', [200, 100, 800, 300]);

%% Cahn-Hilliard dynamics
for istep = 1:nstep
    ttime = ttime + dtime;
    phi1=phi(:,:,1);
    phi2=phi(:,:,2);

    phi1k = fft2(phi1);   
    phi2k = fft2(phi2);   
    phik = cat(3,phi1k,phi2k);

    % Bulk chemical potentials
    mub1 = log(phi1)-log(1-phi1-phi2) + chi11.*phi1 + chi12.*phi2;
    mub1k = fft2(mub1); 
    mub2 = log(phi2)-log(1-phi1-phi2) + chi22.*phi2 + chi12.*phi1;
    mub2k = fft2(mub2); 
        
    % Gradient contributions
    Lphi1k=-k2.*Ldimension^2.*(kappa11.*phi1k + kappa12.*phi2k);
    Lphi1=real(ifft2(Lphi1k));
    Lphi2k=-k2.*Ldimension^2.*(kappa22.*phi2k + kappa12.*phi1k);
    Lphi2=real(ifft2(Lphi2k));

    mu1k=mub1k+Lphi1k;
    mu2k=mub2k+Lphi2k;
    muk=cat(3,mu1k,mu2k);

    % Chemical-potential gradients
    dmudx=real(ifft2(1i*kx.*muk));  
    dmudy=real(ifft2(1i*ky.*muk));
    
    dPix=phi(:,:,1).*dmudx(:,:,1)+phi(:,:,2).*dmudx(:,:,2);
    dPiy=phi(:,:,1).*dmudy(:,:,1)+phi(:,:,2).*dmudy(:,:,2);

    J1x=phi(:,:,1).*(dmudx(:,:,1)-dPix);
    J1y=phi(:,:,1).*(dmudy(:,:,1)-dPiy);

    J2x=phi(:,:,2).*(dmudx(:,:,2)-dPix);
    J2y=phi(:,:,2).*(dmudy(:,:,2)-dPiy);

    dJ1k=1i*kx.*fft2(J1x)+1i*ky.*fft2(J1y);
    dJ2k=1i*kx.*fft2(J2x)+1i*ky.*fft2(J2y);
    dJk=cat(3,dJ1k,dJ2k);

    % Semi-implicit spectral update
    Nonlineark=Ldimension^2.*dJk+A*Ldimension^4.*k4.*phik;
    phik = (Nonlineark.*dtime+phik)./ (1.0 + A.*k4.*Ldimension^4.*dtime);

    phi = real(ifft2(phik));
    phi = cat(3, phi(:,:,1), phi(:,:,2), 1-phi(:,:,1)-phi(:,:,2));

    % Store snapshots and update the figure
    if ttime>(tdraw-dtime/2)
        tdraw=tdraw+savedata_time;
        num=round(ttime/savedata_time);
        num

        phi_store(num,:,:,:)=phi;
        
        % Client concentration map
        subplot(131)
        imagesc(x, y, phi(:,:,2));
        colorbar; 
        shading interp; 
        axis equal tight;
        title(['Client volume fraction \phi_2, t=', num2str(ttime)]);
        
        % Scaffold concentration map
        subplot(132)
        imagesc(x, y, phi(:,:,1));
        colorbar; 
        shading interp;
        axis equal tight;
        title(['Scaffold volume fraction \phi_1, t=', num2str(ttime)])

        % Combined concentration map
        subplot(133)
        image(phi);
        colorbar;
        shading interp; 
        axis equal tight;
        pause(0.05);
        drawnow; 
        
        phi2_new=phi(:,:,2);

        max_delta = max(abs(phi2_new(:) - phi2_old(:)));
        max_delta

        phi2_old=phi(:,:,2);

    end
end

%% Save final state
case_id = 4; % Set to 1, 2, 3, or 4 for the current simulation
case_files = {
    'case1_no_phase_separation.mat'
    'case2_interfacial_adsorption.mat'
    'case3_interfacial_phase_separation.mat'
    'case4_interfacial_phase_separation.mat'
};

save(case_files{case_id}, 'phi');
fprintf('Final state saved to: %s\n', case_files{case_id});