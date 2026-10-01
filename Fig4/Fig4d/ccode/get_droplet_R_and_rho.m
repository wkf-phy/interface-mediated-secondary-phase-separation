function [R, rho, labeled_mask] = get_droplet_R_and_rho(phi)
    phi1 = phi(:,:,1);
    phi2 = phi(:,:,2);
    
    %% Identify droplets
    phi1_min = min(phi1(:));
    phi1_max = max(phi1(:));
    threshold = phi1_min + (phi1_max - phi1_min) / 2.0;
    
    % Label connected regions
    mask = phi1 > threshold;
    [L, ~] = bwlabel(mask, 4);
    
    % Merge droplets connected across the left and right boundaries
    leftLabelsPosition = FindBoundaryMiddleIndex(L(:,1)); 
    for i = leftLabelsPosition' 
        leftlabel = L(i, 1);
        rightlabel = L(i, end);
        if rightlabel ~= 0 && leftlabel ~= 0
            L(L == rightlabel) = leftlabel;
        end
    end
    
    % Merge droplets connected across the top and bottom boundaries
    topLabelsPosition = FindBoundaryMiddleIndex(L(1,:)); 
    for i = topLabelsPosition'
        toplabel = L(1, i);
        bottomlabel = L(end, i);
        if bottomlabel ~= 0 && toplabel ~= 0
            L(L == bottomlabel) = toplabel;
        end
    end
    
    % Relabel droplets consecutively
    unique_labels = unique(L);
    unique_labels(unique_labels == 0) = []; 
    num_features = length(unique_labels);
    labeled_mask = zeros(size(L));
    for k = 1:num_features
        labeled_mask(L == unique_labels(k)) = k;
    end
    
    %% Calculate droplet radii
    R = zeros(num_features, 1);
    props = regionprops(labeled_mask, 'Area');
    for i = 1:num_features
        R(i) = sqrt(props(i).Area / pi);
    end
    
    %% Define droplet influence regions
    all_droplets_mask = (labeled_mask > 0);
    [~, nearest_idx] = bwdist(all_droplets_mask);
    influence_map = labeled_mask(nearest_idx);

    %% Calculate interfacial adsorption
    rho = zeros(num_features, 1);
    for i = 1:num_features
        safe_zone = (influence_map == i);
        phi1_temp=phi1(safe_zone);
        phi1_min = min(phi1_temp(:));
        phi1_max = max(phi1_temp(:));

        % Interfacial weighting function
        W = 4*(phi1 - phi1_min) .* (phi1_max - phi1)./(phi1_max-phi1_min).^2;
        W_local = W(safe_zone);
        phi2_local = phi2(safe_zone);

        rho(i) = sum(W_local .* phi2_local) /2/pi /R(i); % 
    end
end
