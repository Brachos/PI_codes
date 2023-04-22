function [A_h_i] = ellipse_area(Stringer, Fuselage)
    %ELLIPSE AREA computes the swept area between each stringer at each
    %considered cross-sectioN;
    % -----
    %
    % Syntax:
    %   [A_h_i] = ellipse_area(Stringer) computes the swept area [m²] between
    %   each stringer for each considered cross-section of the fuselage.
    %
    % Inputs:
    %   
    %   Stringer:   Stringer.nb: Number of stringers in the cross-section
    %               Stringer.position_y: array of the y position of the
    %               stringer [m]
    %               Stringer.position_z: array of the z position of the
    %               stringer[m]
    %               Stringer.position_x: array of the x position of the
    %               stringer[m]
    %               Stringer.max_load: Maximum load that the stringer can
    %               undergo [N]
    %               Stringer.I_zz: array of the zz-inertia for each
    %               cross-section [m^4]
    %               Stringer.I_yy: array of the yy-inertia for each
    %               cross-section [m^4]
    %               Stringer.I_yz: array of the yz-inertia for each
    %               cross-section [m^4]
    %   Output:     A_h_i: an array containing the swept area between each
    %               stringer for each considered cross-section of the fuselage[m²]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_h_i = zeros(1, Stringer.nb);
theta = 2*pi/Stringer.nb;
    for j=1:Stringer.nb/4
        A_h_i(j) = 0.5*Fuselage.a_max*Fuselage.b_max*(atan((Fuselage.a_max*tan(theta*(j)))/Fuselage.b_max) - atan((Fuselage.a_max*tan(theta*(j-1)))/Fuselage.b_max));    
    end
    A_h_i(Stringer.nb/4 + 1:Stringer.nb/2) = flip(A_h_i(1:Stringer.nb/4));
    A_h_i(Stringer.nb/2 + 1:Stringer.nb) = A_h_i(1:Stringer.nb/2);
end