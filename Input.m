cg = [4.4824 0 0];
% Wing:       
Wing.S = 119.5*0.09290304;      %Wing surface [m²]
%Wing.C_D = 0.023;              %Wing drag coefficient
Wing.MAC = 4.5*0.3048;          %Wing mean aerodynamic chord [m]
%Wing.C_M = ;                   %Wing pitching moment coefficient [Nm]
Wing.l_L = ;                    %Wing lift lever arm [N]
Wing.l_D = ;                    %Wing drag lever arm [m]
Wing.aoa_L = ;                  %Angle of attack between the wing lift lever arm and the chord [rad]
Wing.aoa_D = ;                  %Angle of attack between the wing drag lever arm and the chord [rad]

%Empennage:  
Empennage.S = 29.8*0.09290304;  %Empennage surface [m²]
%Empennage.C_D = ;              %Empennage drag coefficient
Empennage.MAC = 1.1842;         %Empennage mean aerodynamic chord [m]
%Empennage.C_M = ;              %Empennage pitching moment coefficient [Nm]
Empennage.b = 2.17;             %Empennage span [m]
Empennage.l_L = ;               %Empennage lift lever arm [m]
Empennage.T = ;                 %Thrust placed on Empennage [N]
Empennage.l_T = ;               %Empennage thrust lever arm [m]
Empennage.l_D = ;               %Empennage drag lever arm [m]
Empennage.aoa_L = ;             %Angle of attack between the empennage lift lever arm and the chord [rad]
Empennage.aoa_D = ;             %Angle of attack between the empennage drag lever arm and the chord [rad]
Empennage.aoa_T = ;             %Angle of attack between the thrust lift lever arm and the chord [rad]

%Fin:
Fin.AR = 3;                     %Fin aspect ratio [-]
Fin.S = 12.7*0.09290304;        %Fin surface [m²]
Fin.l = ;                       %Fin lever arm [m]

%Flight:
Flight.rho = 0.4594;            %Density [kg/m³]
Flight.V = 211.9616;            %Flight speed [m/s]
Flight.Mach = 0.7;              %Flight Mach number
Flight.aoa = 0;                 %Flight angle of attack [rad]