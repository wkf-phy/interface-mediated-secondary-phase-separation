function con = initializePhi(Nx, Ny, c0)
%INITIALIZEPHI Create a nearly uniform concentration field.
%   Nx, Ny: grid dimensions.
%   c0: target mean concentration.
%   con: concentration field with uniform random perturbations.

    format long; 
    NxNy = Nx * Ny;
    noise = c0/10;  
    
    % Add a random perturbation
    con = zeros(Nx, Ny);
    for i = 1:Nx
        for j = 1:Ny
            con(i, j) = c0 + noise * (0.5 - rand);
        end
    end
end 