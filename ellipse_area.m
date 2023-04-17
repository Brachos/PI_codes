function [A_h_i] = ellipse_area(Stringer)
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
A_h_i = zeros(length(Stringer.positition_y(:,1)), Stringer.nb);
theta = 2*pi/Stringer.nb;
for i=1:length(Stringer.positition_y(:,1))
    for j=1:Stringer.nb/4
        A_h_i(i, j) = 0.5*Stringer.a(i)*Stringer.b(i)*(atan((Stringer.a(i)*tan(theta*(j)))/Stringer.b(i)) - atan((Stringer.a(i)*tan(theta*(j-1)))/Stringer.b(i)));    
    end
    A_h_i(i, Stringer.nb/4 + 1:Stringer.nb/2) = flip(A_h_i(i, 1:Stringer.nb/4));
    A_h_i(i, Stringer.nb/2 + 1:Stringer.nb) = A_h_i(i, 1:Stringer/2);
end
end