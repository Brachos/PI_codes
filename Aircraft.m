% Script that couple all codes together and determine if the aircraft is
% stable or not

%% Parameters
AOA = 0.75*pi/180; %angle of attack
ARw = 7; %Wing Aspect ratio
TAPw = 0.3; %Wing tapper ratio
M = 0.7; %Mach number
Altitude = 30000; %Cruise altitude [feet]

S_F = 0; %Fin area
l_F = 0; %Fin moment arm
le = 0.7; %length of the engine
MT = zeros(3,1); %vector of the total moment 3 different directions for the max weight
Mt = zeros(3,1); %vector of the total moment 3 different directions for the min weight
x_e = 4; %position of the engine inlet
Nelem = 8; % number of differents elements, of different mass
           % (1.Wing;2.Fuselage;3.Tail;4.Engines;5.Landing gears;
           % 6.Payload;7.Fuel;8.Installed weight)
MTOW = 4471; %sum(W); %Maximum Take-Off Weight 

%% Speed
[speed,rho] = speed(Altitude,M);
V_c = speed;

%% Fuselage

D_f_max = 8.169271e-01; %[m]
h_f_max = 2.888274e-01; %[m]
l_f = 5.723540e+00; %[m]

%% Wing

[bw,Sw,CLw_alpha,CLw,CD,D,cw_root,cw_tip,cw_MAC,xw_AC,yw_AC,Vw_fuel,Lambda_LE] = wing(M,Altitude,MTOW);

% bw        = wing span [m]
% Sw        = surface of the wings [m?]
% CLw_alpha = wing 3D lift coefficient slope [1/rad]
% CLw       = wing 3D lift coefficient [-]
% CD        = Drag coefficient [-]
% D         = Drag [N]
% cw_root   = wing root chord [m]
% cw_tip    = wing tip chord [m]
% cw_AC     = wing Mean Aerodynamic Chord (MAC) [m]
% xw_AC     = wing x-coordinate of the Aerodynamic Centre [m]
% yw_AC     = wing y-coordinate of the Aerodynamic Centre [m]
% Vw_fuel   = fuel volume that can be stocked in the wings [l]

xw_cg = 0.4*cw_AC;  %for the wing [35%wMAC;42%wMAC] [m]
%% Fuselage
[D_f_max,a_el,b_el,l_f,V_f]=fuselage_design(MTOW,Vw_fuel);

% a and b are the dimensions of the elliptical cross-section. V_f is the
% volume of the fuselage. 
%% V-Tail
cg_pos = 2;
l_cg = cg_pos;
[S_h,S_v,c_root_h,c_tip_h,c_root_v,c_tip_v, angle, l] = v_tail(MTOW,...
    D_f_max,h_f_max,V_c,cw_MAC,Lambda_LE,Sw, cg_pos,l_f,l_cg,bw);
%%%%%%%%% ENTRY %%%%%%%%%%
% MTOW      = Mass of airplane
% D_f_max   = maximal diameter of fuselage
% h_f_max   = hight of airplaine fuselage (side view)
% V_c       = cruise speed
% cw_MAC    = mean aerodynamic chord
% Lambda_LE = sweep angle of leading edge of the wings
% Sw        = surface of the wings
% cg_pos    = position of center of gravity
% l_f       = fuselage length
% l_cg      = length between center of gravity and nose
% bw        = wing span
%%%%%%%%%% EXIT %%%%%%%%%%%
% S_h       = horizontal surface of the tail
% S_v       = vertical surface of the tail
% c_root_h  = chord at the root of horizontal tail
% c_tip_h   = chord at the tip of horizontal tail
% c_root_v  = chord at the root of vertical tail
% c_tip_v   = chord at the tip of vertical tail
% angle     = dihedral angle of the v-tail
% l         = length between wing ac and tail ac

S_T = S_h; %Tailplane area
l_T = l; %Tail moment arm
hrc = c_root_h; %horizontal tail root chord
htc = c_tip_h; %horizontal tail tip chord 
hTR = c_tip_h/c_root_h; %horizontal tail taper ratio
hMAC = hrc*(2/3)*((1+hTR+hTR^2)/(1+hTR)); %Horizontal Tail Main Aerodynamic Chord
V_hT = S_h*l_T/Sw*hMAC; %Tail volume ratio
xcg_h = 0.3*hMAC; %for the horizontal tail [m]

vrc = c_root_v; %vertical tail root chord
vtc = c_tip_v; %vertical tail tip chord 
vTR = c_tip_v/c_root_v; %vertical tail taper ratio
vMAC = vrc*(2/3)*((1+vTR+vTR^2)/(1+vTR)); %vertical Tail Main Aerodynamic Chord
V_vT = S_v*l_T/Sw*vMAC; %Tail volume ratio
xcg_v = 0.3*vMAC; %for the vertical tail [m]

%% 
%V_F = S_F*l_F/S*c__; %fin volume ratio

%% Weight
[W_wing, W_V, W_fuselage, W_landing_gear, W_installed_weight, W_payload, W_FS, W_tot] = mass(hMAC,vMAC,S_h,S_v,angle,V_hT,V_vT);
W_engine = 140;
W = [W_wing;W_fuselage;W_V;W_engine;W_landing_gear;W_payload;W_FS;W_installed_weight]; %vector of all the different weights (or mass)
                       % (1.Wing;2.Fuselage;3.Tail;4.Engines;5.Landing gears;
                       % 6.Payload;7.Fuel?)
minW = sum(W)-W(6)-W(7); %minimum weight (or minimum mass)
MTOW = sum(W);
%% Center of gravity
x_wv = 4.95; %distance between the wac and the vac
xcg_e = 0.37*le; %for the engine [30%le;45%le] [m]
xcg_f = 0.44*l_f; %for the fuselage [40%L;48%L] [m]
xcg_l= 2; %for the landing gears
xcg_p = 2; %for the payload
x_w = 2; %position of the wings
x_t = 5;
y_wmac = 2; %position of the wing mac along y
y_tmac = 1; %position of the tail mac along y
x_wLE = x_w+sin(Lambda_LE)*y_wmac; %position of the wing leading edge at the mac
x_tLE = x_t+sin(angle)*y_tmac;
thick = 0; %thickness of the plane
zcg_w = 0.07*thick; %for the wing [0.05*thick;0.10*thick] [m]
xarm = [xw_cg+x_wLE;xcg_f;xcg_h+x_tLE;xcg_e+x_e;xcg_l;xcg_p;0;0];
yarm = [0;0;0;0;0;0;0;0]; %symetric
zarm = [zcg_w;0;0;0;0;0;0;0];
% vector of all the different arms corresponding to the different
% weights ; arm = horizontal dimension from the nose of the plane
% to the cg of the mass ; origin : nose
cgT = zeros(3,1); %coordinates of the center of gravity maximum weight
cgt = zeros(3,1); %coordinates of the center of gravity minimum weight
X_cg=0;
for i=1:Nelem
    MT(1) = MT(1) +(W(i)*xarm(i));
    MT(2) = MT(2) +(W(i)*yarm(i));
    MT(3) = MT(3) +(W(i)*zarm(i));
    if(i==6||i==7)
        continue;
    else
        Mt(1) = Mt(1) +(W(i)*xarm(i));
        Mt(2) = Mt(2) +(W(i)*yarm(i));
        Mt(3) = Mt(3) +(W(i)*zarm(i));
    end
end
for i=1:3
    cgT(i) = MT(i)/MTOW; 
    cgt(i) = Mt(i)/minW;
    % Total moment/Total weight ; we can compute a range for the cg with
    % the lower weight and the greater weight different coordinates of the 
    % cg for the max weight
end

%% Neutral point
lwt = 4.95; %Horizontal distance between the wing ac and the tail ac
zwt = 1; %Vertical distance bewteen the wing ac and the tail ac
h0 = 0; %Position of the aerodynamic center
CL = 0;
a = CLw_alpha; %CL_alpha wing
a1 = 0; %CL_alpha tail
r = lwt/(bw/2); 
m = zwt/(bw/2);
de_dAOA = 0.4; %Variation de l'angle epsilon en fonction de l'angle d'attaque
X_np = h0 + V_hT*a1/a*(1-(de_dAOA));

%% Aerodynamic center
D_ac = 0.26*(M-0.4)^2.5; %Delta X_ac ; aerodynamic center
X_c4 = 1/4*cw_AC; %position of the quarter-chord
X_ac = X_c4 + D_ac*sqrt(S);

%% Static margin
k = X_np - X_cg;
%K = -dC_m/dC_Lw;

% -dC_m/dalpha ; static margin is the difference between the position of 
% the neutral point and the position of the cg or the derivative of the 
% coefficient of the total moment ; -C_malpha/C_Lalpha
% 'Certification authorities specify that k >= 0.05