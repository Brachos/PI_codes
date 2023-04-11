% Script that couple all codes together and determine if the aircraft is
% stable or not
clear
%% Figures settings
clc
close all
feature('DefaultCharacterSet','UTF8');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultTextFontsize',13);
set(groot, 'defaultAxesFontsize',13);
set(groot, 'defaultLegendFontsize',13);
set(groot, 'defaultLegendLocation','best');
set(0, 'DefaultLineLineWidth', 1.8);
pt = 12;
%% Parameters
AOA = 2.5; % AOA where the drag is minimum or cl/cd is maximum [deg]
ARw = 7; %Wing Aspect ratio
TAPw = 0.3; %Wing tapper ratio
M = 0.7; %Mach number
Altitude = 30000; %Cruise altitude [feet]
pound = 2.20462262; %kg to lbs
feet = 3.28; %m to ft
inche = 39.37; %m to in
N2Lbf = 0.224809; %N to Lbf factor

MT = zeros(3,1); %vector of the total moment 3 different directions for the max weight
Mt = zeros(3,1); %vector of the total moment 3 different directions for the min weight
Nelem = 10; %number of differents elements, of different mass
% (1.Fuselage;2.Wing;3.Tail;4.Engines+Installed_Weight;5.First Landing
% gears;6.Second Landing Gears;7.Payload;8.Fuel+Installed_Weight;9.System)
MTOW = 4429; %sum(W); [kg] %Maximum Take-Off Weight (Converged, first approx --> 4471)

%% Speed
[speed1] = speed(Altitude,M);
% disp('rho is');
% disp(rho);
rho = 0.48;
V_c = speed1;

%% Wing
[bw,Sw,CLw_alpha,CDw_alpha,CLw,CDw,D,cw_root,cw_tip,cw_MAC,xw_AC,yw_AC,Vw_fuel,sweep,c,alpha_L0,tap,cl_alphaw,theta_tip, A] = wing(M,Altitude,0.95*MTOW,AOA);
% disp('Sw =');
% disp(Sw);
% bw        = wing span [m]
% Sw        = surface of the wings [m?]
% CLw_alpha = wing 3D lift coefficient slope [1/rad]
% CLw       = wing 3D lift coefficient [-]
% CD        = Drag coefficient [-]
% D         = Drag [N]
% cw_root   = wing root chord [m]
% cw_tip    = wing tip chord [m]
% cw_MAC    = wing Mean Aerodynamic Chord (MAC) [m]
% xw_AC     = wing x-coordinate of the Aerodynamic Centre [m]
% yw_AC     = wing y-coordinate of the Aerodynamic Centre [m]
% Vw_fuel   = fuel volume that can be stocked in the wings [l]
% theta_tip = twist angle in RADIANS

xw_cg = 0.4*cw_MAC;  %for the wing [35%wMAC;42%wMAC] [m]

%% Fuselage
[D_f_max,a_el,b_el,l_f,V_f]=fuselage_design(MTOW,Vw_fuel);
% a and b are the dimensions of the elliptical cross-section, semi-axes, a = long axe horizontal, b = small axe vertical. 
% V_f is the volume of the fuselage. 
%h_f_max = 2.888274e-01; %[m]
%l_f = 7; %[m]

%% V-Tail
cg_pos = 4.38;% ? revoir absolument !!!!
l_cg = cg_pos;
[S_tail,Sh_tail,Sv_tail,c_root_tail,c_tip_tail, dihedral_angle, l_arm, CL_tail, Lambda_T,...
    b_tail, bv_tail, bh_tail, W_tail, Cn_beta_Ah, V_vf, hight_root, hight_tip,...
    rudder_chord_root, rudder_chord_tip, rudder_chord, S_rudder, AR_T] = v_tail(MTOW,...
    2*b_el,cw_MAC,sweep*180/pi,Sw,l_f,l_cg,bw);
%%%%%%%%% INPUTS %%%%%%%%%%
% MTOW      = Mass of airplane
% D_f_max   = maximal diameter of fuselage
% h_f_max   = height of airplaine fuselage (side view)
% V_c       = cruise speed
% cw_MAC    = mean aerodynamic chord
% sweep     = sweep angle of leading edge of the wings
% Sw        = surface of the wings
% cg_pos    = position of center of gravity
% l_f       = fuselage length
% l_cg      = length between center of gravity and nose
% bw        = wing span
%%%%%%%%%% OUTPUTS %%%%%%%%%%%
% S_tail            = total surface of the tail
% Sh_tail           = horizontal surface of the tail
% Sv_tail           = vertical surface of the tail
% S_rudder          = Surface of the rudder
% c_root_tail       = chord at the root of the tail
% c_tip_tail        = chord at the tip of the tail
% b_tail            = total span of the tail along itself
% bv_tail           = vertical span of the tail
% bh_tail           = horizontal span of the tail
% dihedral_angle    = dihedral angle of the v-tail in RADIANS
% l_arm             = length between wing ac and tail ac
% Lambda_T          = sweep angle of the tail
% W_tail            = weight of the tail in KILOGRAMS
% hight_root        = hight of the begining of the rudder
% hight_tip         = hight of the end of the rudder
% rudder_chord      = proportion of the chord used for rudder
% rudder_chord_root = rudder chord at the root of the rudder
% rudder_chord_tip  = rudder chord at the tip of the rudder
% AR_T              = aspect ratio of the tail

S_T = Sh_tail; %Tailplane area
l_T = l_arm; %Tail moment arm
hrc = c_root_tail; %horizontal tail root chord
htc = c_tip_tail; %horizontal tail tip chord 
hTR = c_tip_tail/c_root_tail; %horizontal tail taper ratio
hMAC = hrc*(2/3)*((1+hTR+hTR^2)/(1+hTR)); %Horizontal Tail Main Aerodynamic Chord
V_hT = Sh_tail*l_T/(Sw*cw_MAC); %Tail volume ratio
xcg_h = 0.3*hMAC; %for the horizontal tail [m]
xac_h = 0.365*hMAC;

% Vertical
vrc = c_root_tail; %vertical tail root chord
vtc = c_tip_tail; %vertical tail tip chord 
vTR = c_tip_tail/c_root_tail; %vertical tail taper ratio
vMAC = vrc*(2/3)*((1+vTR+vTR^2)/(1+vTR)); %vertical Tail Main Aerodynamic Chord
xcg_v = 0.3*vMAC; %for the vertical tail [m]
xac_v = 0.365*vMAC; 

cl_alphaT = (1.4-1)/(0+4)*180/pi; % Approximated [1/rad]
syms y
c_tail    = (1-2*y/b_tail)*c_root_tail + (2*y/b_tail)*c_tip_tail;
c_MAC_tail = double(2/S_tail*int(c_tail^2,0,b_tail/2));

% fprintf('Yaw coefficient derivative of the Aircraft whithout the fin: %.2dm\n',Cn_beta_Ah);

% Prints
fprintf('Tail surface : %.2dft?\n',S_tail);
fprintf('Tail horizontal span : %.2dm\n',bh_tail);
fprintf('Tail vertical span : %.2dm\n',bv_tail);
fprintf('Dihedral angle (degrees) : %.2d\n',dihedral_angle*180/pi);
fprintf('Surface ratio : %.2d\n',Sh_tail/Sw);
fprintf('Rudder surface : %.2d\n', S_rudder);
fprintf('Total weight of the tail : %.2dkg\n', W_tail);

%% Weight
[W_wing, W_fuselage, W_landing_gear_nose, W_landing_gear_main, W_installed_engine, W_payload, W_FS, W_fuel, W_system, W_tot, W_subsyst, W_sensors] = mass(MTOW,bw,cw_root,cw_tip,l_arm);
W = [W_fuselage;W_wing;W_tail;W_installed_engine;W_landing_gear_nose;W_landing_gear_main;W_payload;W_fuel+W_FS;W_subsyst;W_sensors]; 
%vector of all the different weights (or mass)
% (1.Fuselage;2.Wing;3.Tail;4.Engines+Installed_Weight;5.First Landing
% gears;6.Second Landing Gears;7.Payload;8.Fuel+Fuel system;9.Subsystem;10.Sensors)
minW = sum(W)-W(7)-W(8)+W_FS; %minimum weight (or minimum mass)
MTOW = sum(W);
PayW = sum(W)-W(8)+W_FS;
FW = sum(W)-W(7);

%% Center of gravity
le = 0.7; %length of the engine
x_e = l_f-le; %position of the engine inlet [m]
x_wv = l_arm; %distance between the wac and the vac [m]
xcg_e = 0.37*le; %for the engine [30%le;45%le] [m]
xcg_f = 0.44*l_f; %for the fuselage [40%L;48%L] [m]
xcg_l1= 1.5; %for the first landing gears
xcg_l2= 4.8; %for the second landing gears
xcg_p = 5.2; %for the payload
% xcg_s = 3.6; %for the system (radar...)
xcg_sub = 5; %for the subsystems
xcg_sen = 1.5; %for the sensors
x_w = 2.88; %position of the wings
x_t = l_f-c_root_tail; %position of the tail
xcg_fuel = 3.9; %for the fuel
y_wmac = yw_AC; %position of the wing mac along y
y_tmac = 1; %position of the tail mac along y
syms y
y_hmac = double(2/Sh_tail*int(c*y,0,bh_tail/2));
y_vmac = double(2/Sv_tail*int(c*y,0,bv_tail/2));
x_wLE = x_w+sin(36.3361*pi/180)*y_wmac; %position of the wing leading edge at the mac
x_tLE = x_t+sin(41.3361*pi/180)*y_tmac;
x_hLE = x_t+sin(41.3361*pi/180)*y_hmac;
x_vLE = x_t+sin(41.3361*pi/180)*y_vmac;
thick = 0; %thickness of the plane
zcg_w = 0.07*thick; %for the wing [0.05*thick;0.10*thick] [m]
xarm = [xcg_f;xw_cg+x_wLE;xcg_h+x_tLE;xcg_e+x_e;xcg_l1;xcg_l2;xcg_p;xcg_fuel;xcg_sen;xcg_sub];
yarm = [0;0;0;0;0;0;0;0;0;0]; %symetric
zarm = [zcg_w;0;0;0;0;0;0;0;0;0];
% (1.Fuselage;2.Wing;3.Tail;4.Engines+Installed_Weight;5.First Landing
% gears;6.Second Landing Gears;7.Payload;8.Fuel+Installed_Weight;9.Sensors;10.Subsystems)
% vector of all the different arms corresponding to the different
% weights ; arm = horizontal dimension from the nose of the plane
% to the cg of the mass ; origin : nose
cgT = zeros(3,1); %coordinates of the center of gravity maximum weight
cgt = zeros(3,1); %coordinates of the center of gravity minimum weight
cgp = zeros(3,1); %coordinates of the center of gravity with only payload no fuel
cgf = zeros(3,1); %coordinates of the center of gravity with only fuel no payload
X_cg=0;
MT = zeros(1,3);
Mt = zeros(1,3);
Mp = zeros(1,3);
Mf = zeros(1,3);
for i=1:Nelem
    MT(1) = MT(1) +(W(i)*xarm(i));
    MT(2) = MT(2) +(W(i)*yarm(i));
    MT(3) = MT(3) +(W(i)*zarm(i));
    if(i==7)
        Mp(1) = Mp(1) +(W(i)*xarm(i));
        Mp(2) = Mp(2) +(W(i)*yarm(i));
        Mp(3) = Mp(3) +(W(i)*zarm(i));
        continue;
    elseif(i==8)
        Mf(1) = Mf(1) +(W(i)*xarm(i));
        Mf(2) = Mf(2) +(W(i)*yarm(i));
        Mf(3) = Mf(3) +(W(i)*zarm(i));
        continue;
    else
        Mt(1) = Mt(1) +(W(i)*xarm(i));
        Mt(2) = Mt(2) +(W(i)*yarm(i));
        Mt(3) = Mt(3) +(W(i)*zarm(i));
        Mp(1) = Mp(1) +(W(i)*xarm(i));
        Mp(2) = Mp(2) +(W(i)*yarm(i));
        Mp(3) = Mp(3) +(W(i)*zarm(i));
        Mf(1) = Mf(1) +(W(i)*xarm(i));
        Mf(2) = Mf(2) +(W(i)*yarm(i));
        Mf(3) = Mf(3) +(W(i)*zarm(i));
    end
end
for i=1:3
    cgT(i) = MT(i)/MTOW; 
    cgt(i) = Mt(i)/minW;
    cgp(i) = Mp(i)/PayW;
    cgf(i) = Mf(i)/FW;
    % Total moment/Total weight ; we can compute a range for the cg with
    % the lower weight and the greater weight different coordinates of the 
    % cg for the max weight
end
fprintf('Cg from nose with only fuel, no payload : %.2dm\n',cgf(1));
fprintf('Cg from nose fully loaded : %.2dm\n',cgT(1));
fprintf('Cg from nose empty : %.2dm\n',cgt(1));
fprintf('Cg from nose with only payload, no fuel : %.2dm\n',cgp(1));
h = (cgT(1)-x_wLE)/cw_MAC; % Position of the cg in the case of the MTOW
h2 = (cgt(1)-x_wLE)/cw_MAC; % Position of the cg in the case of the empty aircraft
hp = (cgp(1)-x_wLE)/cw_MAC;
hf = (cgf(1)-x_wLE)/cw_MAC;

%% Lift coefficient
rho_mat = 0.48; %[kg/m^3]
% L = MTOW*9.81;
% CL = L/(Sw*1/2*V_c^2*rho_mat);
% C_L = (CL-CLw)*Sw/S_h;
CL = CLw + CL_tail*Sh_tail/Sw;
L = CL*Sw*1/2*V_c^2*rho_mat;

%% Neutral point
lwt = l_arm; %Horizontal distance between the wing ac and the tail ac
zwt = 1; %Vertical distance bewteen the wing ac and the tail ac
h0 = 0.37; %Position of the aerodynamic center
%CL = 0;
a = CLw_alpha; %CL_alpha wing
a1 = 4; %CL_alpha tail
r = lwt/(bw/2); 
m = zwt/(bw/2);
de_dAOA = 0.45; %Variation de l'angle epsilon en fonction de l'angle d'attaque
de_dAOA1 = 1.75*a/(pi*7*(2*lwt*0.3/bw)^(1/4)*(1+m));
hn = h0 + V_hT*a1/a*(1-(de_dAOA1)); %position of the neutral point in %MAC
x_hn = hn*cw_MAC+x_wLE; %position of the neutral point from nose
fprintf('Neutral point from nose fully loaded : %.2dm\n', x_hn);

%% Aerodynamic center
D_ac = 0.26*(M-0.4)^2.5; %Delta X_ac ; aerodynamic center
X_c4 = 1/4*cw_MAC; %position of the quarter-chord
X_ac = X_c4 + D_ac*sqrt(Sw);
wAC = xw_AC + x_wLE; %Position of the wing aerodynamic center from nose
hAC = xac_h + x_hLE;
vAC = xac_v + x_vLE;
X_w = h*cw_MAC-xw_AC;

%% Pitching moment equation
static_stability = h-(h0+V_hT*a1/a*(1-(de_dAOA1))); %in the case of the MTOW dCm/dalpha
static_stability2 = (h2 - h0)-V_hT*a1/a*(1-(de_dAOA)); %in the case of the empty aircraft

%% Directional stability - Torenbeek
Z_w = b_el; %vertical distance from the wing root quarter chord to the fuselage center line (positive downward)
av_area = pi*a_el*b_el; %average fuselage cross section area
d = sqrt(av_area/0.7854);
% A = 7; %aspect ratio of the wing ? (=3.5 for the tail)
ds_db = -0.276+3.06*Sv_tail/Sw*1/(1+cosd(Lambda_T))+0.4*Z_w/d+0.009*A; % sidewash derivative w.r.t. the yaw angle
V_v = 1; %vertical tailplane airspeed
V = 1; %airspeed
n_v = (V_v/V)^2;
CL_alphaT = 1.1/10 * sin(dihedral_angle);

%% Directional stability (Elsevier)
h_f = 0.5;
Cl_beta_T = - V_vf*h_f/l_arm*CL_alphaT;

%% Derivatives
% Ci-dessous ? revoir !! 
Ixz = 2952;               % Inertia product kg m^2
Ix = 33898;               % Roll moment of inertia kg m^2
Iy = 165669;              % Pitch moment of inertia in kg m^2
Iz = 189496;              % Yaw moment of inertia kg m^2
Inertia = table(Ix, Iy, Iz, Ixz);
V0 = V_c;

% Longitudinal stability
[Long_derivatives,Long_modes] = long_dyn_stab(MTOW,...
    a_el,b_el,bw,Sw,CLw_alpha,rho,V_c,ARw,M,Altitude,CL_alphaT,Sh_tail,...
    de_dAOA1,static_stability,AOA,alpha_L0,l_f,l_cg,sweep,...
    cl_alphaw,l_arm*feet,V_hT,cw_MAC*feet,X_w*feet,cw_root*feet,xw_AC*feet,Inertia);
writetable(Long_derivatives, 'longitudinalStab.txt');

% Lateral stability
[LatDeriv, LatDimDeriv, LatModes, A_lat] = lat_dyn_stab(a_el, b_el, bw, sweep, A, ...
    Sv_tail, Sw, V_vf, dihedral_angle, CLw, l_f, cw_root, V_f, cgT(1), ...
    c_root_tail, bv_tail, Lambda_T, wAC, cw_MAC, theta_tip, M, V_c, AR_T,...
    Sh_tail, c_MAC_tail, CL_tail, bh_tail, x_w, AOA, Inertia, rho, MTOW);
writetable(LatDeriv, 'lateralStab.txt');
% Cn_beta -- Per RADIANS (// Nv in slides)
% Cl_beta -- Per RADIANS (// Lv in slides)
% Clp seems good         (// Lp in slides)

% Criterium for stability, eq. see slide 34 course 02 Flight dynamic :
% x_dot = A*x + B*u, and eigen value of A should have a negative real part
% to ensure stability !
% Ci-dessous ? revoir !! 

%% Static margin
kf = hn - hf;
fprintf('Static margin fuel no payload is about : %.2dm\n',kf);
k = hn - h;
fprintf('Static margin fully loaded is about : %.2dm\n',k);
k2 = hn - h2;
fprintf('Static margin empty is about : %.2dm\n',k2);
kp = hn - hp;
fprintf('Static margin payload no fuel is about : %.2dm\n',kp);
%Range 
x1 = (hn-0.05)*cw_MAC+x_wLE;
x2 = (hn-0.2)*cw_MAC+x_wLE;
%K = -dC_m/dC_Lw;
%K = -static_stability

% -dC_m/dalpha ; static margin is the difference between the position of 
% the neutral point and the position of the cg or the derivative of the 
% coefficient of the total moment ; -C_malpha/C_Lalpha
% 'Certification authorities specify that k >= 0.05

% %% Polar CD_vs_CL
% AOA_vector = -18:0.1:18;
% for i=1:length(AOA_vector)
% [~,~,~,~,CLw,CD,~,~,~,~,~,~,~,~,~] = wing(M,Altitude,MTOW,AOA_vector(i));
% CL_vector(i) = CLw;
% CD_vector(i) = CD;
% CL_CD(i) = CL_vector(i)/CD_vector(i);
% 
% deriv(i) = 0.5/CL_vector(i)*pi*0.8*ARw;
% if abs(deriv(i)*CD_vector(i)-CL_vector(i)) < 0.005
%     CL_opt = CL_vector(i);
%     num = i;
% end
% end
% 
% figure(3)
% plot(CD_vector,CL_vector)
% hold on
% plot(CD_vector,deriv(num)*CD_vector)
% xlim([0 0.2]);
% ylim([0 2]);
% 
% % % Tangente
% % for i = 1 : length(CD_vector)-1
% %     deriv(i) = (CL_vector(i+1)-CL_vector(i))/(CD_vector(i+1)-CD_vector(i));
% %     cd_deriv(i) = (CD_vector(i+1)+CD_vector(i))/2;
% % end
% % 
% % 
% % 
% % figure(3)
% % plot(cd_deriv,deriv)
% 
% figure1 = figure(1);
% clf;
% set(figure1,'defaulttextinterpreter','latex');
% hold on;
% p1=plot(CD_vector,CL_vector,'linewidth', 2, 'MarkerSize', 20');
% xlabel('CD [-]')
% ylabel('CL [-]')
% xlim([0 0.2]);
% ylim([0 2]);
% box on
% figure2 = figure(2);
% clf;
% set(figure2,'defaulttextinterpreter','latex');
% hold on;
% p2=plot(AOA_vector,CL_CD,'linewidth', 2, 'MarkerSize', 20');
% xlabel('AOA [deg]')
% ylabel('CL/CD [-]')
% box on
% % figure;
% % p3=plot(AOA_vector,CD_vector);
% % xlabel('AOA');
% % ylabel('CD');
% figure(3)
% plot(CD_vector,CL_vector)
% hold on
% plot(CD_vector,deriv(num)*CD_vector)
% xlim([0 0.2]);
% ylim([0 2]);
% 
% tail plot
p1 = tand(Lambda_T)*hight_root + rudder_chord_root/rudder_chord - rudder_chord_root;
p2 = tand(Lambda_T)*hight_tip + rudder_chord_tip/rudder_chord - rudder_chord_tip;
figure
plot([0 0 c_root_tail tand(Lambda_T)*bv_tail;...
    c_root_tail tand(Lambda_T)*bv_tail tand(Lambda_T)*bv_tail+c_tip_tail tand(Lambda_T)*bv_tail+c_tip_tail],...
    [0 0 0 bv_tail; 0 bv_tail bv_tail bv_tail],'color',[0 0.4470 0.7410]);
hold on
plot([0 0;-bh_tail/2/3 bh_tail/2/3]-0.2,[0 0; bv_tail/3 bv_tail/3]+0.5,'color',[0.8500 0.3250 0.0980])
hold on 
plot([p1 p1 p2; ...
    p1 + rudder_chord_root p2 p2+rudder_chord_tip],...
    [hight_root hight_root hight_tip; hight_root hight_tip hight_tip],'color',[0.9290 0.6940 0.1250])
title('V-tail geometry')
axis equal
%% 