function[Stringer, Fuselage] = cross_section_geometry(Fuselage, Stringer)
    %CROSS_SECTION_GEOMETRY computes the geometry of each cross-section of the fuselage and places the stringers.
    % -----
    %
    % Syntax:
    %   [] = cross_section_geometry(a_min, a_max, b_min, b_max, x_min, x_max) computes the geometry of each
    %       considered cross-section of the fuselage and the position (x,
    %       y,z) of the strings.
    %
    % Inputs:
    % Fuselage:     Fuselage.x_min: x-position of the minimum cross-section (at tail)[m]
    %               Fuselage.a_min: Minimum semi major axis of the elliptical fuselage (at the tail) [m]
    %               Fuselage.b_min: Minimum semi minor axis of the elliptical fuselage (at the tail) [m]
    %               Fuselage.x_min: x-position of the maximum cross-section (at tail)[m]
    %               Fuselage.a_max: Semi major axis of the fuselage at section A [m]
    %               Fuselage.b_max: Semi minor axis of the fuselage at section A [m]
    %               Fuselage.x_cs: x-position of the considered cross-section [m]
    %               Fuselage.a: Semi major axis of the considered cross-section [m]
    %               Fuselage.b: Semi-minor axis of the considered cross-section [m]
    %               Fuselage.A_h: Area of the cross-sections[m²]
    %               Fuselage.L: Total length of the fuselage [m]
    %               Fuselage.W: Total weight of the fuselage [N]
    %               Fuselage.x: x_positions of the points of the fuselage
    %               cross-sections [m]
    %               Fuselage.y: y_positions of the points of the fuselage
    %               cross-sections [m]
    %               Fuselage.z: z_positions of the points of the fuselage
    %               cross-sections [m]
    % Stringer:     Stringer.nb: number of stringers on the cross-section
    %               Stringer.position_y: y_position of the stringers [m]
    %               Stringer.position_z: z_position of the stringers [m]
    %               Stringer.position_x: x_position of the stringers [m]
    %               Stringer.delta_y: variation of the y_position of the
    %               stringers over the fuselage
    %               Stringer.delta_z: variation of the z_position of the
    %               stringers over the fuselage
    %               Stringer.max_load: tensile Yield Strength [N/m]
    %               Stringer.tau_max: shear strength [N/m]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fuselage.x = Fuselage.x_max*ones(1, 1e4);
Fuselage.y = linspace(-Fuselage.a_max, Fuselage.a_max, 1e4);
Fuselage.z = Fuselage.b_max.*sqrt(1-(Fuselage.y(1, :)./Fuselage.a_max).^2);

%stringers 
Stringer.nb = 24;
theta = 2*pi/Stringer.nb;   %the stringers are uniformly placed

%Benjamin Hardy (2023). ellipse_points 
%(https://www.mathworks.com/matlabcentral/fileexchange/82059-ellipse_points), 
%MATLAB Central File Exchange. Retrieved April 19, 2023. 
[y, z] = ellipse_points(Fuselage.a_max, Fuselage.b_max, Stringer.nb);
Stringer.position_y = y';
Stringer.position_z = z';
[y_min, z_min] = ellipse_points(Fuselage.a_min, Fuselage.b_min, Stringer.nb);
% Stringer.position_y = zeros(1, Stringer.nb); 
% Stringer.position_z = zeros(1, Stringer.nb); 
Stringer.position_x = Fuselage.x_max*ones(1, Stringer.nb);
% y_min= zeros(1, Stringer.nb); 
% z_min = zeros(1, Stringer.nb); 

% for j = 1:Stringer.nb
%     %maximum section
%     r = sqrt(Fuselage.a_max^2*Fuselage.b_max^2/(Fuselage.b_max^2*cos((j-1)*theta)^2 + Fuselage.a_max^2*sin((j-1)*theta)^2));
%     Stringer.position_y(j) = r*cos((j-1)*theta);
%     Stringer.position_z(j) = r*sin((j-1)*theta);
%     %minimum section
%     r = sqrt(Fuselage.a_min^2*Fuselage.b_min^2/(Fuselage.b_min^2*cos((j-1)*theta)^2 + Fuselage.a_min^2*sin((j-1)*theta)^2));
%     y_min(j) = r*cos((j-1)*theta);
%     z_min(j) = r*sin((j-1)*theta);
% end

Stringer.delta_y = (y_min - Stringer.position_y)./(Fuselage.x_min-Fuselage_l_cst - Fuselage.x_max);
Stringer.delta_z = (z_min - Stringer.position_z)./(Fuselage.x_min-Fuselage_l_cst - Fuselage.x_max);
end
