function MiddlePositions=FindBoundaryMiddleIndex(boundary)

    % Identify nonzero labels along the boundary
    lab = unique(boundary);
    boundaryLabels=lab(lab~=0); 
    MiddlePositions = zeros(length(boundaryLabels), 1); 
    
    % Find the midpoint of each labeled boundary segment
    for i = 1:length(boundaryLabels)
        positions = find(boundary == boundaryLabels(i));  
        MiddlePositions(i) = ceil((positions(1) + positions(end)) / 2); 
    end
end