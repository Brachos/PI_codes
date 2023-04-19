function[Stringer] = stringer(Stringer, T_fuselage, M_fuselage)
    %STRINGER computes the minimum area of the stringers for each considered
    %cross-section of the fuselage . It also computes the inertia for each
    %section.
    % -----
    %
    % Syntax:
    %   [Stringer.B] = stringer(Stringer, T_fuselage, M_fuselage) copmutes an
    %   array with the minimum area of the strings [m²] for each cross-section of the fuselage. 
    %   The corresponding inertia are computed as well [m^4].
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
    %   T_fuselage: T_fuselage.X: shear stress along the x-axis for each
    %               considered cross-section and each point of the
    %               enveloppe [N]
    %               T_fuselage.Y: shear stress along the y-axis for each
    %               considered cross-section and each point of the
    %               enveloppe [N]
    %               T_fuselage.Z: shear stress along the z-axis for each
    %               considered cross-section and each point of the
    %               enveloppe [N]
    %   M_fuselage: M_fuselage.X: moment along the x-axis for each
    %               considered cross-section and each point of the
    %               enveloppe [Nm]
    %               M_fuselage.Y: moment along the y-axis for each
    %               considered cross-section and each point of the
    %               enveloppe [Nm]
    %               M_fuselage.Z: moment along the z-axis for each
    %               considered cross-section and each point of the
    %               enveloppe [Nm]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nb_pts_maneuver = length(T_fuselage.X(1, :));
I_zz = sum(Stringer.position_y.^2);
I_yy = sum(Stringer.position_z.^2);
I_yz = sum(Stringer.position_y.*Stringer.position_z);

% minimal boom area for each cross-section
B = zeros(nb_pts_maneuver, Stringer.nb);
for j = 1:nb_pts_maneuver
    sigma_xx_B = @(y, z) ((I_zz*M_fuselage.Y(j) + I_yz*M_fuselage.Z(j)) * z - (I_yz*M_fuselage.Y(j) + I_yy*M_fuselage.Z(j)) * y ) / (I_yy*I_zz - I_yz^2) +  T_fuselage.X(j)/Stringer.nb;
    for t = 1: Stringer.nb
        B(j, t) = abs(sigma_xx_B(Stringer.position_y(t), Stringer.position_z(t)))/Stringer.max_load;
    end
end

Stringer.B = max(B,[],'all');
Stringer.I_yy = Stringer.B*I_yy;
Stringer.I_zz = Stringer.B*I_zz;
Stringer.I_yz = Stringer.B*I_yz;

end
