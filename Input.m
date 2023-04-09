% This script contains all the aircraft data that is needed to compute the
% aerodynamic and structural loads.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cg = [4.4824 0 0];  %Aircraft center of gravity position with respect to the nose

% Aerodynamic center position of the components with respect to the nose
y_ac_wing = 0.5;
y_ac_empennage = 0;
y_ac_fin = 0.7;
y_ac_thrust = 0;

x_ac_wing = 4.1363;
x_ac_empennage = 7.9864;
x_ac_fin = 7.2526;
x_ac_thrust = 7.3540;

%Position of the center of gravity of all the components & their weight
%(1.Fuselage;2.Wing;3.Tail;4.Engines+Installed_Weight;5.First Landing gears;6.Second Landing Gears;7.Payload;8.Fuel+Fuel system;9.Subsystem;10.Sensors)
xarm =[3.6324 4.2821 8.0799 7.8144 1.5000 4.8000 5.0000 4.0000 1.5000 6.0000]; %[m]
W = 1.0e+03 * [0.6743 0.1784 0.0477 0.4079 0.0287 0.0986 0.2050 2.4407 0.1007 0.1034]; %[kg]
    
XX = 1;
%All positions are computed with respect to the nose of the aircraft
%The back of the wing is referred as section A 

%Aircraft:
Aircraft.W = sum(W)*9.18;       %Aircraft weight [N]
Aircraft.I_theta = XX;            %Aircraft inertia [kg.m²]
Aircraft.C_DB = XX;               %Aicraft body drag [N]
Aircraft.l_DB = XX;               %Aicraft drag lever arm [m]

% Wing:       
Wing.S = 119.5*0.09290304;                                  %Wing surface [m²]
Wing.W = 29.8*9.81;                                         %Wing weight [N]
Wing.AR = 7;                                                %Wing aspect ratio [-]
Wing.MAC = 1.07;                                            %Wing mean aerodynamic chord [m]
Wing.l = sqrt((y_ac_wing-cg(2))^2 + (x_ac_wing-cg(1))^2);   %Wing lever arm [m]
Wing.aoa = atan2(y_ac_wing - cg(2), x_ac_wing - cg(1));     %Angle of attack between the wing lever arm and the chord [rad]
Wing.LE = 3.95;                                             %Position of the leading edge [m]
Wing.root_chord = 1.5056;                                   %Root chord [m]
Wing.aoa_fuselage = 1.5 * pi/180;                           %Angle of attack between the wing and the fuselage [rad]


%Empennage:  
Empennage.S = 29.8*0.09290304;                                          %Empennage surface [m²]
Empennage.W = 29.8*9.81;                                                %Empennage weight [N]
Empennage.MAC = 0.75;                                                   %Empennage mean aerodynamic chord [m]
Empennage.b = 2.17;                                                     %Empennage span [m]
Empennage.l = sqrt((y_ac_empennage-cg(2))^2 + (x_ac_empennage-cg(1))^2);%Empennage lever arm [m]
Empennage.T = 8210;                                                     %Thrust placed on Empennage [N]
Empennage.l_T = sqrt((y_ac_thrust - cg(2))^2 + (x_ac_thrust - cg(1))^2);%Empennage thrust lever arm [m]
Empennage.aoa = atan2(y_ac_empennage - cg(2), x_ac_empennage - cg(1));  %Angle of attack between the empennage lever arm and the chord [rad]
Empennage.aoa_T = atan2(y_ac_thrust - cg(2), x_ac_thrust - cg(1));      %Angle of attack between the thrust lift lever arm and the chord [rad]
Empennage.ac = 8.42;                                                    %Aerodynamic centre of the empennage [m]
Empennage.aoa_fuselage = 0;                                             %Angle of attack between the empennage and the fuselage [rad]

%Fin:
Fin.AR = 3;                                                 %Fin aspect ratio [-]
Fin.S = 12.7*0.09290304;                                    %Fin surface [m²]
Fin.l = sqrt((y_ac_fin - cg(2))^2 + (x_ac_fin - cg(1))^2);  %Fin lever arm [m]
Fin.ac = 8.18;                                              %Aerodynamic centre of the fin [m]

%Flight:
Flight.rho = 0.4594;            %Density [kg/m³]

%Engine:
Engine.W = 243*9.81;            %Engine weight [N]
Engine.cg = xarm(4);            %Center of gravity of the engine [m]

%Sensors:
Sensors.W = W(10)*9.81;         %Sensors weight [N]
Sensors.cg = xarm(10);          %Center of gravity of the sensors [m]

%Rear landing gear:
Rear_land_gear.W = W(6)*9.81;   %Rear landing gear weight [N]
Rear_land_gear.cg = xarm(6);    %Center of gravity of the Rear landing gear [m]

%Payload:
Payload.W = W(7)*9.81;          %Payload weight [N]
Payload.cg = xarm(7);           %Center of gravity of the Payload [m]

%Fuselage:  
Fuselage.a_min = XX;            %Minimum semi major axis of the elliptical fuselage (at the tail) [m]
Fuselage.b_min = XX;
Fuselage.a_max = XX;            %Semi major axis of the fuselage at section A [m]
Fuselage.b_max = XX;
Fuselage.L = 8.49;              %Total length of the fuselage [m]
Fuselage.W = W(1)*9.81;         %Total weight of the fuselage [N]
    
%Tail:
Tail.cg = xarm(3);              %Position of the cg of the tail [m]
Tail.W = W(3)*9.81;             %Weight of the v-tail [N]
    
