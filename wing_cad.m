function [coord_cad] = wing_cad(boom_root,boom_tip, stringers_root, stringers_tip, nb_sections,airf_root)

coord_cad.root = zeros(stringers_root.nb,3);
coord_cad.tip = zeros(stringers_root.nb,3);
coord_cad.section1 = zeros(stringers_root.nb,3);
coord_cad.section2 = zeros(stringers_root.nb,3);
coord_cad.section3 = zeros(stringers_root.nb,3);

%for the root
coord_cad.root(:,1) = [stringers_root.XZ_up1(1,:) boom_root.XZ_up(1,1) stringers_root.XZ_up2(1,:) boom_root.XZ_up(1,2) stringers_root.XZ_low1(1,:) boom_root.XZ_low(1,1) stringers_root.XZ_low2(1,:) boom_root.XZ_low(1,2)];
coord_cad.root(:,2) = [stringers_root.XZ_up1(2,:) boom_root.XZ_up(2,1) stringers_root.XZ_up2(2,:) boom_root.XZ_up(2,2) stringers_root.XZ_low1(2,:) boom_root.XZ_low(2,1) stringers_root.XZ_low2(2,:) boom_root.XZ_low(2,2)];

%for the tip
coord_cad.tip(:,1) = [stringers_tip.XZ_up1(1,:) boom_tip.XZ_up(1,1) stringers_tip.XZ_up2(1,:) boom_tip.XZ_up(1,2) stringers_tip.XZ_low1(1,:) boom_tip.XZ_low(1,1) stringers_tip.XZ_low2(1,:) boom_tip.XZ_low(1,2)];
coord_cad.tip(:,2) = [stringers_tip.XZ_up1(2,:) boom_tip.XZ_up(2,1) stringers_tip.XZ_up2(2,:) boom_tip.XZ_up(2,2) stringers_tip.XZ_low1(2,:) boom_tip.XZ_low(2,1) stringers_tip.XZ_low2(2,:) boom_tip.XZ_low(2,2)];
coord_cad.tip(:,3) = coord_cad.tip(:,3) + boom_tip.Y(1);

%interpolation for 3 sub airfoils different sections
% Define section points along the chordwise direction
section_points = linspace(0,1,nb_sections+1);

% Define the z-coordinate of each section
section_z = linspace(boom_root.Y(1), boom_tip.Y(1), nb_sections+2);

% Loop through each section
%for i = 1:nb_sections
    % Calculate interpolation factor for current section
    interp_factor = (section_points(1) + section_points(1+1))/2;
%end
    
    % Interpolate stringer and rib coordinates at current section
    coord_cad.section1(:,1) = (1-interp_factor)*coord_cad.root(:,1) + interp_factor*coord_cad.tip(:,1);
    coord_cad.section1(:,2) = (1-interp_factor)*coord_cad.root(:,2) + interp_factor*coord_cad.tip(:,2);
    %coord_cad.section1(:,3) = (1-interp_factor)*coord_cad.root(:,3) + interp_factor*coord_cad.tip(:,3);
    coord_cad.section1(:,3) = section_z(2) + coord_cad.section1(:,3);
    
    % Interpolate for section 2
    interp_factor = (section_points(2) + section_points(2+1))/2;
    coord_cad.section2(:,1) = (1-interp_factor)*coord_cad.root(:,1) + interp_factor*coord_cad.tip(:,1);
    coord_cad.section2(:,2) = (1-interp_factor)*coord_cad.root(:,2) + interp_factor*coord_cad.tip(:,2);
    coord_cad.section2(:,3) = section_z(3) + coord_cad.section2(:,3);
    
    %coord_cad.section2(:,3) = (1-interp_factor)*coord_cad.root(:,3) + interp_factor*coord_cad.tip(:,3);
    
    % Interpolate for section 3
    interp_factor = (section_points(3) + section_points(3+1))/2;
    coord_cad.section3(:,1) = (1-interp_factor)*coord_cad.root(:,1) + interp_factor*coord_cad.tip(:,1);
    coord_cad.section3(:,2) = (1-interp_factor)*coord_cad.root(:,2) + interp_factor*coord_cad.tip(:,2);
    coord_cad.section3(:,3) = section_z(4) + coord_cad.section3(:,3);

    %coord_cad.section3(:,3) = (1-interp_factor)*coord_cad.root(:,3) + interp_factor*coord_cad.tip(:,3);
    
    xx = zeros(length(airf_root.XZ_up(:,1)),1);
    figure
    
plot3(airf_root.XZ_up(:,1), airf_root.XZ_up(:,2), xx,'linewidth', 2, 'MarkerSize', 11')
hold on
axis equal
plot3(airf_root.XZ_low(:,1),airf_root.XZ_low(:,2),xx,'linewidth', 2, 'MarkerSize', 11')
plot3(coord_cad.root(:,1), coord_cad.root(:,2), coord_cad.root(:,3),'o', 'linewidth', 5,'MarkerSize', 2)
plot3(coord_cad.tip(:,1), coord_cad.tip(:,2), coord_cad.tip(:,3),'o', 'linewidth', 5,'MarkerSize', 2)

plot3(coord_cad.section1(:,1), coord_cad.section1(:,2), coord_cad.section1(:,3),'o', 'linewidth', 5,'MarkerSize', 2)
plot3(coord_cad.section2(:,1), coord_cad.section2(:,2), coord_cad.section2(:,3),'o', 'linewidth', 5,'MarkerSize', 2)
plot3(coord_cad.section3(:,1), coord_cad.section3(:,2), coord_cad.section3(:,3),'o', 'linewidth', 5,'MarkerSize', 2)


coord_cad.root(:, [2, 3]) = coord_cad.root(:, [3, 2]);
coord_cad.tip(:, [2, 3]) = coord_cad.tip(:, [3, 2]);
coord_cad.section1(:, [2, 3]) = coord_cad.section1(:, [3, 2]);
coord_cad.section2(:, [2, 3]) = coord_cad.section2(:, [3, 2]);
coord_cad.section3(:, [2, 3]) = coord_cad.section3(:, [3, 2]);
coord_cad.tip(:, 2) = (-1)*coord_cad.tip(:,2);
coord_cad.section1(:, 2) = (-1)*coord_cad.section1(:,2);
coord_cad.section2(:, 2) = (-1)*coord_cad.section2(:,2);
coord_cad.section3(:, 2) = (-1)*coord_cad.section3(:,2);



writematrix(coord_cad.root*1000,'coord_cad_root.txt');
writematrix(coord_cad.tip*1000,'coord_cad_tip.txt');
writematrix(coord_cad.section1*1000,'coord_cad_section1.txt');
writematrix(coord_cad.section2*1000,'coord_cad_section2.txt');
writematrix(coord_cad.section3*1000,'coord_cad_section3.txt');
%axis equal
end

