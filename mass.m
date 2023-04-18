function [W_wing, W_fuselage, W_landing_gear_nose, W_landing_gear_main, W_installed_weight, W_payload, W_FS, W_mass, W_syst, W_tot, W_subsyst, W_sensors] = mass(MTOW,bw,cw_root,cw_tip,l_arm, net_thrust)
%%Based on chapter 10 of the reference book Aicraft design
%%a sytems engineering approach
Mach = 0.7;
pound = 2.20462262; % kg to lbs
feet = 3.28; % m to ft
inche = 39.37; %m to in
lbf = 0.224809; %N to lbf
knot = 1.94384449; %m/s to knot
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

[~,S_w,~,~,~,~,~,~,~,~,~,~,Vw_fuel,sweep,~,~,tap] = wing(Mach,30000,MTOW,2.5);

W_dg = MTOW*pound - 0.45*W_fuel;
N_z = 1.5*3; %Ultimate load factor = 1.5*limit load factor
S_w = S_w*feet^2; %[ft^2]
AR_w = 7;
tc_root = 0.12; %thickness to chord ratio rekative to the airfoil --> SC-0412
S_csw = 0.3 * S_w; % Control surface area [ft^2] APPROXIMATION
ZFW = MTOW*pound - W_fuel; % Zero fuel weight [lb]

% Raymer General aviation p.606
% W_fw = 0.55*W_fuel;
% q = 216; % [lbf/ft^2] at rho = 0.46 kg/m^3 , speed = 211.96 m^2
% W_wing = 0.036*S_w^0.758*W_fw^0.0035*(AR_w/cos(sweep)^2)^0.6*q^0.006*tap^0.04*(100*tc_root/cos(sweep))^-0.3*(N_z*W_dg)^0.49;
% W_wing = W_wing/pound;

% Raymer Cargo approximation p.604
% W_wing = 0.0051 * (W_dg * N_z)^0.557 * S_w^0.649 * AR_w^0.5 * tc_root^-0.4 * (1 + tap)^0.1 * cos(sweep)^-1 * S_csw^0.1; % [lbs]
% W_wing = W_wing/pound; %[kg]

% Noels approximation
W_wing = 4.22*S_w + 1.642*10^-6 * (N_z*bw^3*sqrt(MTOW*pound*ZFW)*(1+2*tap))/(tc_root*cos(sweep)^2*S_w*(1+tap));
W_wing = W_wing/pound;

%% Installed engine weight Raymer
N_E = 1; %[-] nbre of engines
K_E = 2.575; %[N] using metric units --> engine weight factor
W_E = 243*pound; %[kg] weight of each engine

W_installed_weight = K_E*N_E*W_E^0.922; %[lbs]
W_installed_weight = W_installed_weight/pound; %[kg]

%% Fuselage Raymer
[D_f_max,~,~,l_f,~] = fuselage_design(MTOW,Vw_fuel,net_thrust);
L_f = l_f*feet; % [m] fuselage length
D_f = D_f_max*feet; %[ft] fuselage max diameter of the eq. circ. cross-sect??
A_top = L_f*D_f; % Approximation of the fuselage as a cylinder
A_side = A_top;
Swet = 3.4*((A_top+A_side)/2); % Fuselage wetted area [ft^2] Raymer p.235
K_door = 1.12; % =1.12 if aft clamshell door (door for the payload)
K_Lg = 1.12; % =1.12 if fuselage mounted main landing gear
L = L_f; % Fuselage structural length [ft]
lambda = tap; % Wing taper ratio
Kws = 0.75*((1+2*lambda)/(1+lambda)) * (bw/L) * tan(sweep); % Wing sweep factor
D = D_f; % Fuselage structural depth [ft]

% Raymer general aviation approximation
% Lt = l_arm*feet;
% W_press = 0;
% W_fuselage = (0.052 * (Swet)^1.086 * (N_z*W_dg)^0.177 * (Lt^(-0.051)) * (L_f/D_f)^-0.072 * q^0.241 + W_press)/pound; % General aviation Raymer p.606

% Noels approximation
% Dpmax = 73773.90*lbf/feet^2; % Cabin pressure at 2600 m [Pa]-->[lb/ft^2] Conceptual design slide 71
% Nlim = 3; % Limit load factor
% I_p = 1.5*10^-3 * Dpmax * D_f; % Pressure index
% ZFW = MTOW*pound - W_fuel; % Zero fuel weight [lb]
% Ww = W_wing*pound; % Wing weight
% We = W_installed_weight*pound; % Engine weight
% I_b = 1.91*10^-4 * Nlim * (ZFW - Ww - We) * L_f/D_f^2; % Bending index
% if I_p>I_b
%     I_fus = I_p;
% elseif I_p<I_b
%     I_fus = (I_p^2 + I_b^2)/(2*I_b);
% end
% W_fuselage = (1.051 + 0.102*I_fus) * Swet;
% W_fuselage = W_fuselage/pound;

% Raymer Cargo approximation
W_fuselage = 0.3280 * K_door * K_Lg * (W_dg * N_z)^0.5 * L^0.25 * Swet^0.302 * (1+Kws)^0.04 * (L/D)^0.10; % Cargo Raymer p.604 [lbs]
W_fuselage = W_fuselage/pound; %[kg]

%% Nose landing gear Raymer
N_l = 1.5; % nombre de landing gear * 1.5
W_l = MTOW*pound; %in [lbs]
L_n = 0.8*inche; %in [in] hauteaur
Knp = 1; % =1.15 for kneeling gear, 1 otherwise
Nnw = 1; % Number of nosewheels

%Raymer General aviation approx
% W_landing_gear = 0.125*(N_l*W_l)^0.566*(L_n/12)^0.845; %Raymer p.606
% W_landing_gear_nose  = W_landing_gear/pound; %in [kg]

%Raymer Cargo approx p.605
W_landing_gear_nose = 0.032 * Knp * W_l^0.646 * N_l^0.2 * L_n^0.5 * Nnw^0.45;
W_landing_gear_nose = W_landing_gear_nose/pound;

%% Main landing gear Raymer
L_m = 0.8*inche; %in [in]
Kmp = 1; % 1.126 for kneeling gear, 1 otherwise
Nmw = 2; % Number of main wheels
Nmss = 2; % Number of main gear shock struts
Vstall = 242.29*knot; % Stall speed depending on rho,CL,MTOW and Sw [m/s]-->[kt]

% Rayper general aviation approx
% W_landing_gear = 0.095*(N_l*W_l)^0.768*(L_m/12)^0.409;
% W_landing_gear_main = W_landing_gear/pound; %in [kg]

% Raymer cargo approx
W_landing_gear_main = 0.0106 * Kmp * W_l^0.888 * N_l^0.25 * L_m^0.4 * Nmw^0.321 * Nmss^-0.5 * Vstall^0.1;
W_landing_gear_main = W_landing_gear_main/pound;

%% Payload
W_payload = 110; %[kg]

%% System
IESuP = 222/pound; % Initial Estimated Subsystem Payload
IESeP = 27; % Initial Estimated Sensors Payload 
W_syst = IESuP+IESeP;
W_subsyst = IESuP;
W_sensors = IESeP;
%% Total mass
W_tot = W_landing_gear_nose + W_landing_gear_main + W_fuselage + W_FS + W_payload + W_installed_weight + W_wing + W_mass; %add the fuel weight
end