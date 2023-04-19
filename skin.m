function[thickness_min] = skin(Stringer, T_fuselage, M_fuselage, Fuselage)
    %SKIN computes the minimum thickness of the fuselage skin.
    % -----
    %
    % Syntax:
    %   [thickness] = skin(Stringer, T_fuselage, M_fuselage) computes the 
    %   minimum thickness of the skin of the fuselage [m].
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
    %   Outputs:
    %   thickness: the minimum thickness of the skin of the fuselage [m].
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Stress in the web
nb_pts_maneuver = length(T_fuselage.X(1, :));
thickness = zeros(1, nb_pts_maneuver);
A_h_i = ellipse_area(Stringer, Fuselage);

for t=1:nb_pts_maneuver
    P_x = zeros(1, Stringer.nb);
    P_y = zeros(1, Stringer.nb);
    P_z = zeros(1, Stringer.nb);
    for j=1:Stringer.nb
        P_x(j) = ((Stringer.I_zz*M_fuselage.Y(t) + Stringer.I_yz*M_fuselage.Z(t))*Stringer.position_z(j) - (Stringer.I_yz*M_fuselage.Y(t) + Stringer.I_yy*M_fuselage.Z(t))*Stringer.position_y(j)) / (Stringer.I_yy*Stringer.I_zz - Stringer.I_yz^2)*Stringer.B;
        P_y(j) = P_x(j)*Stringer.delta_y(j);
        P_z(j) = P_x(j)*Stringer.delta_z(j);
    end
    T_webb.y(t) = T_fuselage.Y(t) - sum(P_y);
    T_webb.z(t) = T_fuselage.Z(t) - sum(P_z);
    
    %% Shear flux
    q = zeros(1, Stringer.nb); %open shear flux (cut between 1&2)[N/m]
    
    for j=2:Stringer.nb
        q(j) = q(j-1) - (Stringer.I_zz*T_webb.z(t) - Stringer.I_yz*T_webb.y(t))/(Stringer.I_yy*Stringer.I_zz - Stringer.I_yz^2) * Stringer.position_z(j)*Stringer.B;
        q(j) = q(j) - (Stringer.I_yy*T_webb.y(t) - Stringer.I_yz*T_webb.z(t))/(Stringer.I_yy*Stringer.I_zz - Stringer.I_yz^2) * Stringer.position_y(j)*Stringer.B;
    end
    q0 = (-sum(P_z.*Stringer.position_y) + sum(P_y.*Stringer.position_z) - 2*sum(q.*A_h_i))/(2*Fuselage.A_h);   %closed shear flux[N/m]
    q_t = M_fuselage.X(t)/(2*Fuselage.A_h);                                                                 %Twist part of the shear flux[N/m]
    
    q_tot = q + q0*ones(1,Stringer.nb) + q_t*ones(1,Stringer.nb);   %total shear flux[N/m]
    thickness(t) = max(abs(q_tot))/Stringer.tau_max;                %minimum thickness [m]
end
thickness_min = max(thickness);
end
