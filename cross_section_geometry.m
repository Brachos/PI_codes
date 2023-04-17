function[Stringer] = cross_section_geometry(Fuselage, Stringer)
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

%3 considered cross-sections
Fuselage.x = Fuselage.x_cs'.*ones(1, 1e4);
Fuselage.y = zeros(length(Fuselage.a), 1e4);
Fuselage.z = zeros(length(Fuselage.a), 1e4);
for i = 1:length(Fuselage.a)
    Fuselage.y(i, :) = linspace(-Fuselage.a(i), Fuselage.a(i), 1e4);
    Fuselage.z(i, :) = Fuselage.b(i).*sqrt(1-(Fuselage.y(i, :)./Fuselage.a(i)).^2);
end

%stringers 
Stringer.nb = 24;
theta = 2*pi/Stringer.nb;   %the stringers are uniformly placed

Stringer.position_y = zeros(length(Fuselage.x_cs), Stringer.nb); 
Stringer.position_z = zeros(length(Fuselage.x_cs), Stringer.nb); 
Stringer.position_x = Fuselage.x_cs'.*ones(1, Stringer.nb);

for i = 1:length(Fuselage.a)
    for j = 1:Stringer.nb   
        r = sqrt(Fuselage.a(i)^2*Fuselage.b(i)^2/(Fuselage.b(i)^2*cos((j-1)*theta)^2 + Fuselage.a(i)^2*sin((j-1)*theta)^2));
        Stringer.position_y(i, j) = r*cos((j-1)*theta);
        Stringer.position_z(i, j) = r*sin((j-1)*theta);
    end
end

Stringer.delta_y = (Stringer.position_y(end, :) - Stringer.position_y(1, :))./(Stringer.position_x(end, :) - Stringer.position_x(1, :));
Stringer.delta_z = (Stringer.position_z(end, :) - Stringer.position_z(1, :))./(Stringer.position_x(end, :) - Stringer.position_x(1, :));
end
