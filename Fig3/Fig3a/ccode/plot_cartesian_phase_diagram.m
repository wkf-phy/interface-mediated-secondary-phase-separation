function plot_cartesian_phase_diagram(chit11, chit12, chit22, c1, c2)
% Plot a ternary phase diagram in Cartesian coordinates (phi1, phi2).
% chit11, chit12, chit22: reduced interaction parameters.
% c1: optional fixed-phi1 line; use [] to omit.
% c2: retained for compatibility; not used in this function.

    % Convert to pairwise interaction parameters
    chi13 = -chit11 / 2;
    chi23 = -chit22 / 2;
    chi12 = chit12 + chi13 + chi23;
    
    % Composition grid and Flory-Huggins free energy
    n_points = 1000; 
    [phi2, phi3] = meshgrid(linspace(0,1, n_points), linspace(0,1, n_points));
    phi1 = 1 - phi2 - phi3;
    valid_idx = (phi1 >=0) & (phi2 >=0) & (phi3 >=0);
    
    free_energy = zeros(size(phi1));
    for i = 1:numel(phi1)
        if valid_idx(i)     
            phi1(i) = max(phi1(i), 1e-10);
            phi2(i) = max(phi2(i), 1e-10);
            phi3(i) = max(phi3(i), 1e-10);
            free_energy(i) = phi1(i)*log(phi1(i)) + phi2(i)*log(phi2(i)) + phi3(i)*log(phi3(i)) + ...
                             chi12*phi1(i)*phi2(i) + chi13*phi1(i)*phi3(i) + chi23*phi2(i)*phi3(i);
        else    
            free_energy(i) = NaN; 
        end
    end
    
    x_coord = phi1;
    y_coord = phi2;
    
    x_valid = x_coord(valid_idx);
    y_valid = y_coord(valid_idx);
    z_valid = free_energy(valid_idx);
    points = [x_valid, y_valid, z_valid]; 
    
    % Construct the convex hull and select faces by their normal
    K = convhulln(points);
    
    lower_face = [];
    for i = 1:size(K, 1)
        v12 = points(K(i, 2), :) -  points(K(i, 1), :);
        v13 = points(K(i, 3), :) -  points(K(i, 1), :);
        normal = cross(v12, v13);
        normal = normal / norm(normal); 
        if normal(3) > 0.1
            lower_face = [lower_face; K(i, :)];
        end
    end
    
    % Pure-component vertices
    A = [1, 0];   
    B = [0, 1];   
    C = [0, 0];   
    
    % Classify hull faces using a numerical edge-length threshold
    num_faces = size(lower_face, 1);
    all_edge_lengths = zeros(num_faces, 3);
    for i = 1:num_faces
        idx = lower_face(i, :);
        p1 = points(idx(1), 1:2);
        p2 = points(idx(2), 1:2);
        p3 = points(idx(3), 1:2);
        all_edge_lengths(i, 1) = norm(p1 - p2);
        all_edge_lengths(i, 2) = norm(p2 - p3);
        all_edge_lengths(i, 3) = norm(p3 - p1);
    end

    max_len = max(all_edge_lengths(:));
    threshold = max_len / 5; 
    count_long_edges = sum(all_edge_lengths > threshold, 2);
    idx_1phase = find(count_long_edges == 0); 
    idx_2phase = find(count_long_edges == 2); 
    idx_3phase = find(count_long_edges == 3); 
    
    hold on; axis equal;  axis on;
 
    xlim([0, 0.3]);
    ylim([0, 0.2]);
    xlabel('\phi_1', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('\phi_2', 'FontSize', 12, 'FontWeight', 'bold');

    V = points(:, 1:2); 
    
    % Shade the one-, two-, and three-phase regions
    if ~isempty(idx_1phase)
        patch('Faces', lower_face(idx_1phase, :), 'Vertices', V, 'FaceColor', [0.95 0.75 0.75], 'EdgeColor', 'none');
    end
    if ~isempty(idx_2phase)
        patch('Faces', lower_face(idx_2phase, :), 'Vertices', V, 'FaceColor', [0.85 0.93 1.00], 'EdgeColor', 'none');
    end
    if ~isempty(idx_3phase)
        patch('Faces', lower_face(idx_3phase, :), 'Vertices', V, 'FaceColor', [1.00 0.92 0.80], 'EdgeColor', 'none');
    end
    
    % Compute the eigenvalue of the free-energy Hessian
    s_phi2 = phi2; s_phi3 = phi3; s_phi1 = 1 - s_phi2 - s_phi3;
    min_val = 1e-6;
    s_phi1(s_phi1 < min_val) = min_val; 
    s_phi2(s_phi2 < min_val) = min_val; 
    s_phi3(s_phi3 < min_val) = min_val;
    
    H22 = (1 ./ s_phi1) + (1 ./ s_phi2) - 2 * chi12;
    H33 = (1 ./ s_phi1) + (1 ./ s_phi3) - 2 * chi13;
    H23 = (1 ./ s_phi1) + chi23 - chi12 - chi13;

    Tr = H22 + H33;                     
    Det = H22 .* H33 - H23.^2;          
    Discriminant = Tr.^2 - 4 * Det;     
    lambda1 = (Tr + sqrt(Discriminant)) / 2;
    lambda2 = (Tr - sqrt(Discriminant)) / 2;
    min_eig = min(lambda1, lambda2);
    
    invalid_mask = (phi1 + phi2 + phi3 - 1 > 1e-5) | (phi1 < 0) | (phi2 < 0) | (phi3 < 0);
    min_eig(invalid_mask) = NaN;
    
    % Plot the spinodal
    s_x = phi1; s_y = phi2;
    
    [~, ~] = contour(s_x, s_y, min_eig, [0, 0], 'LineColor', '#006400', 'LineWidth', 2, 'LineStyle', '--');
    
    Z_hatch = s_x + s_y; 
    Z_hatch(min_eig > 0 | isnan(min_eig)) = NaN; 
    [~, ~] = contour(s_x, s_y, Z_hatch, 70, 'LineColor', '#008000', 'LineWidth', 0.5, 'LineStyle', '-');
    
    % Optionally mark a fixed-phi1 line and its three-phase boundary intersections
    if ~isempty(c1)
        p_start_x = c1;
        p_start_y = 0;
        p_end_x = c1;
        p_end_y = 1 - c1;

        plot([p_start_x, p_end_x], [p_start_y, p_end_y], 'k--', 'LineWidth', 1.5);
        text(p_start_x, p_start_y - 0.05, ['\phi_1=', num2str(c1)], 'FontSize', 10, 'HorizontalAlignment', 'center');

        if ~isempty(idx_3phase)
            tp_indices = lower_face(idx_3phase, :);
            ThreePhase_Vertices = points(tp_indices(:), 1:2);
            k_3p = convhull(ThreePhase_Vertices(:,1), ThreePhase_Vertices(:,2));
            bound_pts = ThreePhase_Vertices(k_3p, 1:2);
            
            P1 = [p_start_x, p_start_y];
            P2 = [p_end_x, p_end_y];
            
            intersect_pts = [];
            intersect_phi2 = [];
            
            for i = 1:(length(k_3p)-1)
                Q1 = bound_pts(i, :);
                Q2 = bound_pts(i+1, :);
                
                A_mat = [(P2 - P1)', (Q1 - Q2)'];
                if abs(det(A_mat)) > 1e-8 
                    tu = A_mat \ (Q1 - P1)';
                    t = tu(1); u = tu(2);
                    if t >= -1e-5 && t <= 1+1e-5 && u >= -1e-5 && u <= 1+1e-5
                        int_pt = P1 + t * (P2 - P1);
                        intersect_pts = [intersect_pts; int_pt];
                        int_phi2 = int_pt(2);
                        intersect_phi2 = [intersect_phi2; int_phi2];
                    end
                end
            end
        
            % Remove duplicate intersections and label their phi2 values
            if ~isempty(intersect_pts)
                rounded_pts = round(intersect_pts, 4);
                [~, unique_idx] = unique(rounded_pts, 'rows');
                intersect_pts = intersect_pts(unique_idx, :);
                intersect_phi2 = intersect_phi2(unique_idx);
                
                for i = 1:size(intersect_pts, 1)
                    ix = intersect_pts(i, 1);
                    iy = intersect_pts(i, 2);
                    iphi2 = intersect_phi2(i);
                    
                    scatter(ix, iy, 40, 'rs', 'filled', 'MarkerEdgeColor', 'k');
                    text(ix + 0.02, iy, sprintf('\\phi_2=%.3f', iphi2), ...
                        'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');
                end
            end
        end
    end
    
end