function[] = stringer(Stringer, T_fuselage, M_fuselage)
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

I_zz = zeros(1, length(Stringer.positition_y(:,1)));
I_yy = zeros(1, length(Stringer.positition_y(:,1)));
I_yz = zeros(1, length(Stringer.positition_y(:,1)));
for i = 1:length(Stringer.positition_y(:,1))
    I_zz(i) = sum(Stringer.positition_y(i,:).^2);
    I_yy(i) = sum(Stringer.positition_z(i,:).^2);
    I_yz(i) = sum(Stringer.positition_y(i,:).*Stringer.positition_z(i,:)); 
end

% minimal boom area for each point of the load enveloppe and each cross-section
Stringer.B = zeros(1, length(Stringer.positition_y(:,1))); 
Stringer.I_yy = zeros(1, length(Stringer.positition_y(:,1)));
Stringer.I_zz = zeros(1, length(Stringer.positition_y(:,1)));
Stringer.I_yz = zeros(1, length(Stringer.positition_y(:,1)));

for i = 1:length(Stringer.positition_y(:,1))
    B = zeros(length(Stringer.positition_y(:,1)), String.nb);
    for j = 1:length(Stringer.positition_y(:,1)) 
    sigma_xx_B = @(y, z) ((I_zz(i)*M_fuselage.Y(i,j) + I_yz(i)*M_fuselage.Z(i,j)) * z - (I_yz(i)*M_fuselage.Y(i,j) + I_yy(i)*M_fuselage.Z(i,j)) * y ) / (I_yy(i)*I_zz(i) - I_yz(i)^2) +  T_fuselage.X(i,j)/String.nb;    
        for t = 1: String.nb
            B(j, t) = asb(sigma_xx_B(String.position_y(i), String.position_z(i)))/String.max_load;
        end
    end
    Stringer.B(i) = max(B,[],'all');
    Stringer.I_yy(i) = Stringer.B(i)*I_yy;
    Stringer.I_zz(i) = Stringer.B(i)*I_zz;
    Stringer.I_yz(i) = Stringer.B(i)*I_yz;
end
end
