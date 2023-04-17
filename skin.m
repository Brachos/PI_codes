function[thickness_min] = skin(Stringer, T_fuselage, M_fuselage)
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
nb_cross_sections = length(Stringer.positition_y(:,1))-1;
nb_pts_maneuver = length(T_fuselage.X(1, :));
thickness = zeros(1, nb_pts_maneuver);
A_h_i = ellipse_area(Stringer);

for t=1:nb_pts_maneuver
    P_x = zeros(nb_cross_sections, Stringer.nb);
    P_y = zeros(nb_cross_sections, Stringer.nb);
    P_z = zeros(nb_cross_sections, Stringer.nb);
    for i=1:length(nb_cross_sections)
        for j=1:Stringer.nb
            P_x(i,j) = ((Stringer.I_zz(i)*M_fuselage.Y(i, t) + Stringer.I_yz(i)*M_fuselage.Z(i, t))*Stringer.position_z(i, j) - (Stringer.I_yz(i)*M_fuselage.Y(i, t) + Stringer.I_yy(i)*M_fuselage.Z(i, t))*Stringer.postition_y(i,j)) / (Stringer.I_yy(i)*Stringer.I_zz(i) - Stringer.I_yz(i)^2)*Stringer.B(i);
            P_y(i,j) = P_x(i,j)*Stringer.delta_y(j);
            P_z(i,j) = P_x(i,j)*Stringer.delta_z(j);
        end
        T_webb.y(i, t) = T_fuselage.Y(i, t) - sum(P_y(i, :));
        T_webb.z(i, t) = T_fuselage.Z(i, t) - sum(P_z(i, :));
    end
    
    %% Shear flux
    q = zeros(nb_cross_sections, Stringer.nb); %open shear flux (cut between 1&2)[N/m]
    q0 = zeros(nb_cross_sections, 1);          %closed shear flux[N/m]
    q_t = zeros(nb_cross_sections,1);          %Twist part of the shear flux[N/m]
    
    for i=1:nb_cross_sections
        for j=2:Stringer.nb
            q(i, j) = q(i, j-1) - (Stringer.I_zz(i)*T_webb.z(i, t) - Stringer.I_yz(i)*T_webb.y(i, t))/(Stringer.I_yy(i)*Stringer.I_zz(i) - Stringer.I_yz(i)^2) * Stringer.position_z(i, j)*Stringer.B;
            q(i, j) = q(i, j) - (Stringer.I_yy(i)*T_webb.y(i, t) - Stringer.I_yz(i)*T_webb.z(i, t))/(Stringer.I_yy(i)*Stringer.I_zz(i) - Stringer.I_yz(i)^2) * Stringer.position_y(i, j)*Stringer.B;
        end
        q0(i) = (-sum(P_z(i, :).*String.position_y(i, :)) + sum(P_y(i, :).*String.position_z(i, :)) - 2*sum(q(i,:).*A_h_i(i, :)))/(2*Stringer.A_h(i));
        q_t(i) = M_fuselage.X(i, t)/(2*Stringer.A_h(i));
        
    end
    q_tot = q + q0*ones(1,Stringer.nb) + q_t*ones(1,Stringer.nb);   %total shear flux[N/m]
    thickness(t) = max(abs(q_tot))/Stringer.tau_max;                %minimum thickness [m]
end
thickness_min = max(thickness);
end
