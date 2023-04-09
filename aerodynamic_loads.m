function [L_W, L_E, M_fus, F_fin, D_B, M_tail] = aerodynamic_loads(Aircraft, Wing, Empennage, Fin, Flight)
    %AERODYNAMIC_LOADS computes the aerodynamic loads of the aircraft.
    % -----
    %
    % Syntax:
    %   [L_W, L_E, M_fus, F_fin,D_B] = aerodynamic_loads(AOA, n, Aircraft, Wing, Empennage, Fin, Flight) returns an array with
    %   the wing load [N], the empennage load [N], the fuselage moment [Nm], the fin load [N] and the body drag[N].
    %
    % Inputs:
    %   
    %   Aircraf:    Aircraft.W: weight of the aircraft [N]
    %               Aircraft.I_theta: inertia of aircraft around CG
    %               Aircraft.C_DB: body drag coefficient [kg/m²]
    %               Aircraft.lB: body drag lever arm [m]
    %   Wing:       Wing.S: Wing surface [m²]
    %               Wing.W: Wing weight [N]
    %               Wing.AR: Wing aspect ratio [-]
    %               Wing.MAC: Wing mean aerodynamic chord [m]
    %               Wing.l: Wing lever arm [m]
    %               Wing.aoa: Angle of attack between the wing lever
    %               arm and the chord [rad]
    %               Wing.LE: Position of the leading edge [m]
    %               Wing.root_chord: Wiing root chord [m]
    %               Wing.aoa_fuselage: Angle of attack between the wing and
    %               the fuselage center line [rad]
    %   Empennage:  Empennage.S: Empennage surface [m²]
    %               Empennage.W: Empennage weight [N]
    %               Empennage.MAC: Empennage mean aerodynamic chord [m]
    %               Empennage.b: Empennage span [m]
    %               Empennage.l: Empennage lever arm [m]
    %               Empennage.T: Thrust placed on Empennage [N]
    %               Empennage.l_T: Empennage thrust lever arm [m]
    %               Empennage.aoa: Angle of attack between the empennage lever
    %               arm and the chord [rad]
    %               Empennage.aoa_T: Angle of attack between the thrust lift lever
    %               arm and the chord [rad]
    %               Empennage.ac: Aerodynamic centre of the empennage [m]
    %               Empennage.aoa_fuselage: Angle of attack between the empennage and
    %               the fuselage center line [rad]
    %   Fin:        Fin.AR: Fin aspect ratio [-]
    %               Fin.S: Fin surface [m²]
    %               Fin.l: Fin lever arm [m]
    %               Fin.ac: Aerodynamic centre of the fin [m]
    %   Flight:     Flight.rho: Density [kg/m³]
    %               Flight.V: Flight speed [m/s]
    %               Flight.Mach: Flight Mach number 
    %               Flight.aoa: Flight angle of attack [rad]
    %               Flight.n: load factor [m/s²]
    %
    % Outputs:
    %   L_W: Wing load [N]
    %   L_E: Empennage load [N]
    %   M_fus: Fuselage moment [Nm]
    %   F_fin: Fin load [N]
    %   D_B: Body drag [N]
    %   M_tail: Tail moement [Nm]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %%FAR
 psi = deg2rad(15);         %[rad]
 theta_ddot = deg2rad(60);  %[rad/s²]
 
 %%Fin load
 A = 2*Fin.AR;
 a1 = (5.5*A)/(A+2);
 F_fin = 0.5 * Flight.rho * Flight.V^2 * Fin.S * a1 * psi;
 
 %%Drag
 D_W = 0.5 * Wing.C_D * Flight.rho * Flight.V^2 * Wing.S;
 D_E = 0.5 * Empennage.C_D * Flight.rho * Flight.V^2 * Empennage.S;
 D_B = 0.5 * Aircraft.C_DB * Flight.rho * Flight.V^2 * (Wing.S + Empennage.S);
 
 %%Pitching moments
 M_W = Wing.C_M * 0.5 * Flight.rho * Flight.V^2 * Wing.S * Wing.MAC;
 M_E = Empennage.C_M * 0.5 * Flight.rho * Flight.V^2 * Empennage.S * Empennage.MAC;
 
 %%Tailplane torque
 M_tail = 0.00125/(sqrt(1-Flight.Mach^2)) * Flight.rho * Flight.V^2 * Empennage.S * Empennage.b * psi;
 
 %%Total torque
 M_fus = M_tail + F_fin*Fin.l;

 %%Distance with respect of the AOA
 Wing.l_L = Wing.l * cos(Wing.aoa - Flight.aoa);
 Empennage.l_L = Empennage.l * cos(Empennage.aoa - Flight.aoa + Wing.aoa_fuselage - Empennage.aoa_fuselage);
 Aircraft.l_DB = Aircraft.l_DB * cos(Flight.aoa);
 Empennage.l_D = Empennage.l * sin(Empennage.aoa - Flight.aoa + Wing.aoa_fuselage - Empennage.aoa_fuselage);
 Wing.l_D = Wing.l * sin(Wing.aoa - Flight.aoa);
 Empennage.l_T = Empennage.l_T * sin(Empennage.aoa_T - Flight.aoa + Wing.aoa_fuselage - Empennage.aoa_fuselage);
 
 %%System of equations for L_E and L_W
 A = [1 1; Wing.l -Empennage.l];
 b = [Flight.n*Aircraft.W - Empennage.T*sin(Flight.aoa - Wing.aoa_fuselage + Empennage.aoa_fuselage); Aircraft.I_theta*theta_ddot - M_E - M_W - D_B*Aircraft.l_DB + D_E*Empennage.l - D_W*Wing.l + Empennage.T*Empennage.l_T];

 %%Vector containing L_W and L_C
 x = A\b; 
 L_W = x(1);
 L_E = x(2);

end