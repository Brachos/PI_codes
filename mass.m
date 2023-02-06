function [W_wing, W_V, W_fuselage, W_landing_gear, W_installed_weight, W_payload, W_FS, W_mass, W_tot] = mass(hMAC,vMAC,S_h,S_v,angle,V_hT,V_vT,MTOW,bw,cw_root,cw_tip,bh,l)
%%Based on chapter 10 of the reference book Aicraft design
%%a sytems engineering approach
Mach = 0.7;
%% Fuel system weight for transport and fighter aircraft 
%equipped with integral fuel tanks
[W_mass] = fuel_weight();
W_fuel = W_mass*2.20462262; %[lbs] total fuel weight
rho_f = 5.87; %[lb/gal] fuel density
N_E = 1; % nbre of engines
N_t = 3; % 2 wings + tank inside fuselage = nbre of separated fuel tanks

W_FS = 15*N_t^0.5 * (W_fuel/rho_f)^0.333 + 80*(N_E + N_t -1);
W_FS = W_FS/2.20462262; %[kg]

%% Wing 
% assumption : Remote-controlled model
g = 9.81; %[m/s^2] % gravity
n_max = 1.75; % 7.5 if fighter % ultimate load factor
n_ult = 1.5*n_max; % [-]
% S_w = 12.1; % Wing planform area
% MAC = 1.44; % Wing mean aerodynamic chord
% max_t_c_ratio = 14/100; %Maximum thickness-to-chord ratio
% rho_mat = 2711; % Density of construction material, here aerospace Al
% K_rho = 0.00125; % 0.005 if fighter, few light stores under wing
% AR = 7; % aspect ratio
Lambda = 36.34; % quarter chord sweep angle
lambda = 0.3; % taper ratio
% 
% W_wing = S_w*MAC*max_t_c_ratio*rho_mat*K_rho*((AR*n_ult)/cos(Lambda))^0.6 * lambda^0.04 * g;
%% Wing Raymer
S_w = 12.1*3.28^2; % Wing planform area [ft^2]
AR = 7; % aspect ratio
K_dw = 1;
K_vs = 1;
W_dg = MTOW*2.20462262-0.45*W_fuel; %W_fuel in pounds
Nz = 1.5*3; %take-off
Scsw = 1/4*(0.3*cw_root+0.3*cw_tip)*bw*3.28^2; %[ft^2]
W_wing = 0.0103*K_dw*K_vs*(W_dg*Nz)^0.5*S_w^0.622*AR^0.785*0.14*(1+lambda)^0.05*cos(Lambda*pi/180)^-1*Scsw^0.04;
W_wing = W_wing/2.20462262;
%0.14 = t/c_root

%% V-tail CHECK ALL THE VALUES
%assumption : mass V_tail = mass Horizontal-tail + mass Vertical-tail
%H-tail
%K_rho_HT = 0.075; % [-] 0.07 is supersonic fighter %HT density factor
rho_mat = 2711; % Density of construction material, here aerospace Al
K_rho_HT = 0.0175;
S_HT = S_h; % HT exposed planform area
MAC_HT = hMAC; % HT mean aerodynamic chord
t_c_ratio_HT = 0.365; %HT max thickness-to-chord ratio
AR_HT = 4.7; %HT aspect ratio
Lambda_HT = angle; %HT quarter chord sweep angle
lambda_HT = 0.5; %HT taper ratio
V_H = V_hT; %HT volume ratio
elevator_tail_HT = 0.3; %elevator-to-tail chord ratio [0.2;0.4]

W_HT = S_HT*MAC_HT*t_c_ratio_HT*rho_mat*K_rho_HT*(AR_HT/cos(Lambda_HT))^0.6 * lambda_HT^0.04 * V_H^0.3*elevator_tail_HT^0.4*g;

%V-tail
%K_rho_VT = 0.052; %[-] 0.135 if supersonic fighter %VT density factor
K_rho_VT = 0.052;
S_VT = S_v; % VT exposed planform area
MAC_VT = vMAC; % VT mean aerodynamic chord
t_c_ratio_VT = t_c_ratio_HT; %VT max thickness-to-chord ratio
AR_VT = AR_HT; %VT aspect ratio
Lambda_VT = angle; %VT quarter chord sweep angle
lambda_VT = 0.3; %VT taper ratio
V_V = V_vT; %VT volume ratio
rudder_tail_VT = 0.3; %rudder-to-tail chord ratio [0.2;0.4]

W_VT = S_VT*MAC_VT*t_c_ratio_VT*rho_mat*K_rho_VT*(AR_VT/cos(Lambda_VT))^0.6 * lambda_VT^0.04 * V_V^0.2*rudder_tail_VT^0.4*g;

W_V = W_VT + W_HT; %total mass of the V-tail.

%% Horizontal tail Raymer
Fw = 3.28; %[ft]
W_HT = 3.316*(1+Fw/bh)^-2*(W_dg*Nz/1000)^0.26*S_HT^0.806; %[pounds]
W_HT = W_HT/2.20462262;

%% Vertical tail Raymer
Krht = 1;
HT_HV = 1;
Lt = l;
Sr = 0;
W_VT = 0.452*Krht*(1+HT_HV)^0.5*(W_dg*Nz)^0.488*S_VT^0.718*Mach^0.341*Lt^-1*(1+Sr/S_VT)^0.348*AR_VT^0.223*(1+lambda_VT)^0.25*(cos(Lambda_VT))^-0.323;
W_VT = W_VT/2.20462262;

%% Total tail weight

%% Fuselage 
L_f = 5.72; % [m] fuselage length
D_f = 0.82; %[m] fuselage max diameter of the eq. circ. cross-sect??
K_rhof = 0.002; % 0.0075 if fighter fuselage density factor
K_inlet = 1; %inlets not on fuselage

W_fuselage =  L_f * D_f^2 * rho_mat * K_rhof * n_ult^0.25 * K_inlet * g;
%% Landing gear
b = 9.2; %wing span
K_ret = 1.07; %retractable landing gear
K_L = 1; % not navy aircraft
K_LG = 0.435; %0.335 fighter
H_LG = 0.5; %landing gear height % TO BE MODIFIED
n_ult_land = 1;%landing ultimate factor % TO BE MODIFIED
W_L = 4471; % Landing weight %TO BE MODIFIED

W_landing_gear = K_L*K_ret*K_LG*W_L*(H_LG/b)*n_ult_land^0.2;
%% Installed engine weight
N_E = 1; %[-] nbre of engines
K_E = 3; %[N] using metric units
W_E = 140; %[kg] weight of each engine

W_installed_weight = K_E*N_E*W_E^0.9;

%% Payload
W_payload = 205.0237512;

%% Total mass
W_tot = W_landing_gear + W_fuselage + W_FS + W_payload + W_installed_weight + W_V + W_wing + W_mass; %add the fuel weight
end