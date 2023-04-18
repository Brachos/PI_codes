function [W_wing, W_fuselage, W_landing_gear_nose, W_landing_gear_main, W_installed_weight, W_payload, W_FS, W_mass, W_syst, W_tot, W_subsyst, W_sensors] = mass(MTOW,bw,cw_root,cw_tip,l_arm, net_thrust)
%%Based on chapter 10 of the reference book Aicraft design
%%a sytems engineering approach
Mach = 0.7;
pound = 2.20462262; % kg to lbs
feet = 3.28; % m to ft
inche = 39.37; %m to in
% rho_mat = 2711; % Density of construction material, here aerospace Al


%% Fuel system weight for transport and fighter aircraft 
%equipped with integral fuel tanks
[W_mass] = fuel_weight(net_thrust);
W_fuel = W_mass*pound; %[lbs] total fuel weight
rho_f = 6.67632; %[lb/gal] fuel density for a fuel of 800 kg/m^3
N_E = 1; % nbre of engines
N_t = 3; % 2 wings + tank inside fuselage = nbre of separated fuel tanks
Vt = W_fuel/rho_f;

W_FS = 2.49*Vt^(0.726)*(1/(1+0.99))^(0.363)*N_t^(0.242)*N_E^(0.157); %Raymer p.606
W_FS = W_FS/pound; %[kg]

%% Wing (Raymer formula) (british units)

[~,S_w,~,~,~,~,~,~,~,~,~,~,Vw_fuel,sweep] = wing(Mach,30000,MTOW,2.5);
g = 9.81; %[m/s^2] % gravity
W_dg = MTOW*pound - 0.45*W_fuel;
N_z = 1.5*3; %Ultimate load factor = 1.5*limit load factor
S_w = S_w*feet^2; %[ft^2]
AR_w = 7;
tc_root = 0.14; %thickness to chord ratio
lambda_w = 0.3; % Tapper ratio
Lambda_w = sweep; % Sweep angle

% General aviation 
W_fw = 0.55*W_fuel;
q = 10685*pound/feet^2; %[kg/m^2]
W_wing = 0.036*S_w^0.758*W_fw^0.0035*(AR_w/cos(Lambda_w)^2)^0.6*q^0.006*lambda_w^0.04*(100*tc_root/cos(Lambda_w))^-0.3*(N_z*W_dg)^0.49;
W_wing = W_wing/pound;

%% Fuselage
[D_f_max,~,~,l_f,~] = fuselage_design(MTOW,Vw_fuel,net_thrust);
L_f = l_f*feet; % [m] fuselage length
D_f = D_f_max*0.8*feet; %[ft] fuselage max diameter of the eq. circ. cross-sect??
Lt = l_arm*feet;
W_press = 0;

W_fuselage = (0.052*(pi*D_f*L_f+pi*D_f^2/2)^1.086*(N_z*W_dg)^0.177*Lt^-0.051*(L_f/D_f)^-0.072*q^0.241+W_press)/pound; %Raymer p.606
% W_fuselage = 460;
%% Nose landing gear Raymer
N_l = 1.5; % nombre de landing gear * 1.5
W_l = MTOW*pound; %in [lbs]
L_n = 0.8*inche; %in [in] hauteaur

W_landing_gear = 0.125*(N_l*W_l)^0.566*(L_n/12)^0.845; %Raymer p.606
W_landing_gear_nose  = W_landing_gear/pound; %in [kg]

%% Main landing gear Raymer
L_m = 0.8*inche; %in [in]

W_landing_gear = 0.095*(N_l*W_l)^0.768*(L_m/12)^0.409;
W_landing_gear_main = W_landing_gear/pound; %in [kg]

%% Installed engine weight Raymer
N_E = 1; %[-] nbre of engines
K_E = 2.575; %[N] using metric units --> engine weight factor
W_E = 243*pound; %[kg] weight of each engine

W_installed_weight = K_E*N_E*W_E^0.922;
W_installed_weight = W_installed_weight/pound; %[kg]
%% Payload
W_payload = 110; %[kg]

%% System
IESuP = 222/pound; % Initial Estimated Subsystem Payload
IESeP = 228/pound; % Initial Estimated Sensors Payload
W_syst = IESuP+IESeP;
W_subsyst = IESuP;
W_sensors = IESeP;
%% Total mass
W_tot = W_landing_gear + W_fuselage + W_FS + W_payload + W_installed_weight + W_wing + W_mass; %add the fuel weight
end