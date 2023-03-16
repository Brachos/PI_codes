function [L_W, L_E, M_fus, F_fin, D_B] = aerodynamic_loads(AOA, n, W, I_theta, C_DB, l_DB, Wing, Empennage, Fin, Flight)
    %AERODYNAMIC_LOADS returns the equivalent cruise and dive speed at alt.
    %It also plots the placard diagram of the airplane.
    % -----
    %
    % Syntax:
    %   [L_W, L_E, M_fus, F_fin,D_B] = aerodynamic_loads(AOA, n, W, I_theta, C_DB, l_DB, Wing, Empennage, Fin, Flight) returns an array with
    %   the wing load [N], the empennage load [N], the fuselage moment [Nm], the fin load [N] and the body drag[N].
    %
    % Inputs:
    %   AOA: angle of attack [rad]
    %   n: load factor [m/s²]
    %   W: weight of the aircraft [kg]
    %   I_theta: inertia of aircraft around CG
    %   C_DB: body drag coefficient [kg/m²]
    %   l_DB: body drag lever arm [m]
    %   Wing:       Wing.S: Wing surface [m²]
    %               Wing.C_D: Wing drag coefficient
    %               Wing.c: Wing mean aerodynamic chord [m]
    %               Wing.C_M: Wing pitching moment coefficient [Nm]
    %               Wing.l_L: Wing lift lever arm [N]
    %               Wing.l_D: Wing drag lever arm [m]
    %               Wing.aoa_L: Angle of attack between the wing lift lever
    %               arm and the chord [rad]
    %               Wing.aoa_D: Angle of attack between the wing drag lever
    %               arm and the chord [rad]
    %   Empennage:  Empennage.S: Empennage surface [m²]
    %               Empennage.C_D: Empennage drag coefficient
    %               Empennage.c: Empennage mean aerodynamic chord [m]
    %               Empennage.C_M: Empennage pitching moment coefficient
    %               [Nm]
    %               Empennage.b: Empennage span [m]
    %               Empennage.l_L: Empennage lift lever arm [m]
    %               Empennage.T: Thrust placed on Empennage [N]
    %               Empennage.l_T: Empennage thrust lever arm [m]
    %               Empennage.l_D: Empennage drag lever arm [m]
    %               Empennage.aoa_L: Angle of attack between the empennage lift lever
    %               arm and the chord [rad]
    %               Empennage.aoa_D: Angle of attack between the empennage drag lever
    %               arm and the chord [rad]
    %               Empennage.aoa_T: Angle of attack between the thrust lift lever
    %               arm and the chord [rad]
    %   Fin:        Fin.AR: Fin aspect ratio [-]
    %               Fin.S: Fin surface [m²]
    %               Fin.l: Fin lever arm [m]
    %   Flight:     Flight.rho: Density [kg/m³]
    %               Flight.V: Flight speed [m/s]
    %               Flight.Mach: Flight Mach number 
    %               Flight.aoa: Flight angle of attack [rad]
    %
    % Outputs:
    %   L_W: Wing load [N]
    %   L_E: Empennage load [N]
    %   M_fus: Fuselage moment [Nm]
    %   F_fin: Fin load [N]
    %   D_B: Body drag [N]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %%FAR
 psi = deg2rad(15);         %[rad]
 theta_ddot = deg2rad(60);  %[rad/s²]
 
 %%Fin load
 A = 2*Fin.AR;
 a1 = (5.5*A)/(A+2);
 F_fin = 0.5 * Flight_rho * Flight.V^2 * Fin.S * a1 * psi;
 
 %%Drag
 D_W = 0.5 * Wing.C_D * Flight.rho * Flight.V^2 * Wing.S;
 D_E = 0.5 * Empennage*C_D * Flight.rho * Flight.V^2 * Empennage.S;
 D_B = 0.5 * C_DB * Flight.rho * Flight.V^2 * (Wind.S + Empennage.S);
 
 %%Pitching moments
 M_W = Wing.C_M * 0.5 * Flight.rho * Flight.V^2 * Wing.S * Wing.c;
 M_E = Empennage.C_E * 0.5 * Flight.rho * Flight.V^2 * Empennage.S * Empennage.c;
 
 %%Tailplane torque
 M_tail = 0.00125/(sqrt(1-Flight.Mach^2)) * Flight.rho * Flight.V^2 * Empennage.S * Empennage.b * psi;
 
 %%Total torque
 M_fus = M_tail + F_fin*Fin.l;

 %%Distance with respect of the AOA
 Wing.l_L = Wing.l_L * cos(Wing.aoa_L - AOA);
 Empennage.l_L = Empennage.l_L * cos(Empennage.aoa_L - AOA);
 l_DB = l_DB * cos(AOA);
 Empennage.l_D = Empennage.l_D * cos(Empennage.aoa_D - AOA);
 Wing.l_D = Wing.l_D * cos(Wing.aoa_D - AOA);
 Empennage.l_T = Empennage.l_T * cos(Empennage.aoa_T - AOA);
 
 %%System of equations for L_E and L_W
 A = [1 1; Wing.l_L -Empennage.l_L];
 b = [n*W - Empennage.T*sin(AOA); I_theta*theta_ddot - M_E - M_W - D_B*l_DB + D_E*Empennage.l_D - D_W*Wing.l_D + Empennage.T*Empennage.l_T];

 %%Vector containing L_W and L_C
 x = A\b; 
 L_W = x(1);
 L_E = x(2);

end