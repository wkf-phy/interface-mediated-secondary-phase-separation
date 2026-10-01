load("critical_profile.mat")

%% Parameters
c1_scan_list = 0.07:0.01:0.20; 
color_map = lines(length(c1_scan_list));

%% Plot critical scaffold profiles
figure(1);hold on; box on;
for i = 1:length(c1_scan_list)
    plot(critical_profiles(i).r, critical_profiles(i).p1, ...
        'LineWidth', 2, 'Color', color_map(i,:), ...
        'DisplayName', sprintf('\\phi_{1,avg} = %.2f', critical_profiles(i).c1_base));
end

%% Plot critical client profiles
figure(2);hold on; box on;
for i = 1:length(c1_scan_list)
    plot(critical_profiles(i).r, critical_profiles(i).p2, ...
            'LineWidth', 2, 'Color', color_map(i,:), ...
            'DisplayName', sprintf('\\phi_{1,avg} = %.2f', critical_profiles(i).c1_base));
end

%% Calculate critical interfacial density
for i = 1:length(c1_scan_list)
    r=critical_profiles(i).r;
    phi1=critical_profiles(i).p1;
    phi2=critical_profiles(i).p2;

    phi1max=max(phi1);
    phi1min=min(phi1);

    phi1_target = (max(phi1) + min(phi1)) / 2;
    [~, idx] = min(abs(phi1 - phi1_target));
    r_interface = r(idx);
    
    % Calculate the interfacial weighting function and density
    W=4*(phi1-phi1min).*(phi1max-phi1)./(phi1max-phi1min).^2;
    rho=2*pi*trapz(r, W.*phi2.*r)/2/pi/r_interface;
    [r_interface,rho]

    r_l(i)=r_interface;
    rho_l(i)=rho;
end

%% Plot size dependence of critical density
figure(3);
hold on;

line_color = [0.8500, 0.3250, 0.0980]; 

plot(r_l, rho_l, '-o', ...
    'LineWidth', 2, ...
    'Color', line_color, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', line_color, ...
    'MarkerEdgeColor', 'w'); 

xlabel('$r$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Surface density $\rho$', 'Interpreter', 'latex', 'FontSize', 14);

% Axes formatting
ylim([0.01, 0.015]);
set(gca, 'FontSize', 17, 'LineWidth', 1.2,  ...
    'GridColor', [0.5 0.5 0.5]); 

xticks(0.12:0.04:0.28);
yticks(0.01:0.001:0.015);
title('Size Dependence of Critical Density')

%% Export
% exportgraphics(gca, 'FigureS2.pdf', ...
%     'ContentType', 'vector');