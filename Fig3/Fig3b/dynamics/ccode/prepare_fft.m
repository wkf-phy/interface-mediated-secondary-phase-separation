function [kx, ky, k2, k4] = prepare_fft(Nx, Ny, dx, dy)
%PREPARE_FFT Construct Fourier wave numbers and their squared magnitudes.

    format long

    Nx2 = Nx / 2 + 1;
    Ny2 = Ny / 2 + 1;
    
    % Wave-number spacing.
    delkx = (2.0 * pi) / (Nx * dx);
    delky = (2.0 * pi) / (Ny * dy);

    % Wave numbers in FFT ordering.
    for i = 1:Nx2
        fk1 = (i - 1) * delkx;
        kx(i) = fk1;
        kx(Nx - i + 2) = -fk1;
    end

    for j = 1:Ny2
        fk2 = (j - 1) * delky;
        ky(j) = fk2;
        ky(Ny - j + 2) = -fk2;
    end

    % Squared magnitude of each wave vector.
    for i = 1:Nx
        for j = 1:Ny
            k2(i, j) = kx(i)^2 + ky(j)^2;
        end
    end

    k4 = k2.^2;
end