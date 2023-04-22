% This script contains all the aircraft data that is needed to compute the
% aerodynamic and structural loads.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cg = [4.3537 0 0.6];  %Aircraft center of gravity position with respect to the nose

% Aerodynamic center position of the components with respect to the nose
y_ac_wing = 1;
y_ac_empennage = 0.67;
y_ac_fin = 0.8796;
y_ac_thrust = 0.33956;

x_ac_wing = 2.7635;
x_ac_empennage = 7.89;
x_ac_fin = 7.89;
x_ac_thrust = 8.1506;

%Position of the center of gravity of all the components & their weight
%(1.Fuselage;2.Wing;3.Tail;4.Engines+Installed_Weight;5.First Landing gears;6.Second Landing Gears;7.Payload;8.Fuel+Fuel system;9.Subsystem;10.Sensors)
xarm =[3.4358 4.3555 8.1224 7.5942 1.3000 4.8000 4.5000 4.2780 1.800 0.8000]; %[m]
    
W = [445.324212885535 114.260212744617 60.9058010298828 242.800000000000 29.2282889652049 100.518080161912 150 2322.03916274679 100.697506224444 40];%kg   
%All positions are computed with respect to the nose of the aircraft
%The back of the wing is referred as section A 

%Aircraft:
Aircraft.W = sum(W)*9.18;                   %Aircraft weight [N]
Aircraft.I_theta = 13395974401.883000000e-6;%Aircraft inertia [kg.m²]
Aircraft.C_DB = 0.038299150427751;          %Aicraft body drag coefficient[-]
Aircraft.l_DB = 0.1;                        %Aicraft drag lever arm [m] !!!

% Wing:       
Wing.S = 7.986662323478662;                                 %Wing surface [m²]
Wing.W = W(2);                                              %Wing weight [N]
Wing.AR = 7;                                                %Wing aspect ratio [-]
Wing.MAC = 1.171387198829498;                               %Wing mean aerodynamic chord [m]
Wing.l = sqrt((y_ac_wing-cg(3))^2 + (x_ac_wing-cg(1))^2);   %Wing lever arm [m]
Wing.aoa = atan2(y_ac_wing - cg(3), x_ac_wing - cg(1));     %Angle of attack between the wing lever arm and the chord [rad]
Wing.LE = 2.459;                                            %Position of the leading edge [m]
Wing.root_chord = 1.643312976775196;                        %Root chord [m]
Wing.aoa_fuselage = 1.5 * pi/180;                           %Angle of attack between the wing and the fuselage [rad] !!!


%Empennage:  
Empennage.S = 1.502493233419805;                                        %Empennage surface [m²]
Empennage.W = 24.6516*9.81;                                             %Empennage weight [N]
Empennage.MAC = 0.875342478575083;                                      %Empennage mean aerodynamic chord [m]
Empennage.b = 2.398819063189448;                                        %Empennage span [m]
Empennage.l = sqrt((y_ac_empennage-cg(3))^2 + (x_ac_empennage-cg(1))^2);%Empennage lever arm [m]
Empennage.T = 3298.22248722500;                                         %Thrust placed on Empennage [N] 
Empennage.l_T = sqrt((y_ac_thrust - cg(3))^2 + (x_ac_thrust - cg(1))^2);%Empennage thrust lever arm [m]
Empennage.aoa = atan2(y_ac_empennage - cg(3), x_ac_empennage - cg(1));  %Angle of attack between the empennage lever arm and the chord [rad]
Empennage.aoa_T = atan2(y_ac_thrust - cg(3), x_ac_thrust - cg(1));      %Angle of attack between the thrust lift lever arm and the chord [rad]
Empennage.ac = x_ac_empennage;                                          %Aerodynamic centre of the empennage [m]
Empennage.aoa_fuselage = 0;                                             %Angle of attack between the empennage and the fuselage [rad]
Empennage.AR = 4;                                                       %Aspect ratio

%Fin:
Fin.AR = 4;                                                 %Fin aspect ratio [-]
Fin.S = 1.438583224480276;                                  %Fin surface [m²]
Fin.l = sqrt((y_ac_fin - cg(3))^2 + (x_ac_fin - cg(1))^2);  %Fin lever arm [m]
Fin.ac = x_ac_fin;                                              %Aerodynamic centre of the fin [m]

%Flight:
Flight.rho = 0.4594;            %Density [kg/m³]

%Engine:
Engine.W = W(4)*9.81;           %Engine weight [N]
Engine.cg = xarm(4);            %Center of gravity of the engine [m]
Engine.l = 1.3200852;           %length of the engine [m]
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
Fuselage.x_min = 8.1506;        %x-position of the minimum cross-section [m]
Fuselage.a_min = 0.33956;       %Minimum semi major axis of the elliptical fuselage (at the tail) [m]
Fuselage.b_min = 0.33956;       %Minimum semi minor axis of the elliptical fuselage (at the tail) [m]
Fuselage.x_max = 3.959;         %x-position of the cross-section A [m]
Fuselage.a_max = 0.7461;        %Semi major axis of the fuselage at section A [m]
Fuselage.b_max = 0.514;         %Semi minor axis of the fuselage at section A [m]
Fuselage.A_h = pi*Fuselage.a_max*Fuselage.b_max;%Area of the cross-sections[m²]
Fuselage.L = 8.589668283077916; %Total length of the fuselage [m]
Fuselage.W = W(1)*9.81;         %Total weight of the fuselage [N]
    
%Tail:
Tail.cg = xarm(3);              %Position of the cg of the tail [m]
Tail.W = W(3)*9.81;             %Weight of the v-tail [N]
    
