%% Parameters
AOA = 0; %angle of attack
M = 0.7; %Mach number
S_F = 0; %Fin area
l_F = 0; %Fin moment arm
le = 0.7; %length of the engine
L = 5.72; %length of the fuselage
MT = zeros(3,1); %vector of the total moment 3 different directions for the max weight
Mt = zeros(3,1); %vector of the total moment 3 different directions for the min weight
x_e = 0; %position of the engine inlet
Nelem = 0; % number of differents elements, of different mass
           % (1.Wing;2.Fuselage;3.Tail;4.Engines;5.Landing gears;
           % 6.Payload;7.Fuel?)
           
%% Wing
wCL = 0.4; %CL de la wing
b = 9.2; %wingspan [m]
S = 12.1; %wing area S=15 [m]
c_ = S/b; %standard mean chord (smc) [m]
AR = b^2/S; %Aspect Ratio (from 1.5 to 18)
x_wLE = 2; %position of the leading edge of the wMAC
wrc = 2.02; %wing root chord [m]
wtc = 0.61; %wing tip chord [m]
wTR = wtc/wrc; %wing taper ratio
wMAC = wrc*(2/3)*((1+wTR+wTR^2)/(1+wTR)); %Wing Main Aerodynamic Chord (MAC)
xcg_w = 0.4*wMAC;  %for the wing [35%wMAC;42%wMAC] [m]

%% V-Tail
S_T = 0; %Tailplane area
l_T = 0; %Tail moment arm
x_wv = 4.95; %distance between the wac and the vac
vrc = 0; %v-tail root chord
vtc = 0; %v-tail tip chord 
vTR = 0; %v-tail taper ratio
vMAC = vrc*(2/3)*((1+vTR+vTR^2)/(1+vTR)); %V-Tail Main Aerodynamic Chord
V_T = S_T*l_T/S*vMAC; %Tail volume ratio
xcg_v = 0.3*vMAC; %for the v-tail [m]

%% 
%V_F = S_F*l_F/S*c__; %fin volume ratio

%% Weight
W = [0;0;0;140;0;0;0]; %vector of all the different weights (or mass)
                       % (1.Wing;2.Fuselage;3.Tail;4.Engines;5.Landing gears;
                       % 6.Payload;7.Fuel?)
minW = sum(W)-W(6)-W(7); %minimum weight (or minimum mass)
MTOW = 5000; %sum(W); %Maximum Take-Off Weight

%% Center of gravity
xcg_e = 0.37*le; %for the engine [30%le;45%le] [m]
xcg_f = 0.44*L; %for the fuselage [40%L;48%L] [m]
xcg_l= 2; %for the landing gears
xcg_p = 2; %for the payload
thick = 0; %thickness of the plane
zcg_w = 0.07*thick; %for the wing [0.05*thick;0.10*thick] [m]
xarm = [xcg_w+x_wLE;xcg_f;xcg_w+x_wLE+x_wv;xcg_e+x_e;xcg_l;xcg_p;0];
yarm = [0;0;0;0;0;0;0]; %symetric
zarm = [zcg_w;0;0;0;0;0;0];
% vector of all the different arms corresponding to the different
% weights ; arm = horizontal dimension from the nose of the plane
% to the cg of the mass ; origin : nose
cgT = zeros(3,1); %coordinates of the center of gravity maximum weight
cgt = zeros(3,1); %coordinates of the center of gravity minimum weight
X_cg=0;
for i=1:Nelem
    MT(1) = MT(1) +(W(i)*xarm(i));
    MT(2) = MT(2) +(W(i)*xarm(i));
    MT(3) = MT(3) +(W(i)*xarm(i));
    if(i==6||i==7)
        continue;
    else
        Mt(1) = Mt(1) +(W(i)*xarm(i));
        MT(2) = MT(2) +(W(i)*xarm(i));
        MT(3) = MT(3) +(W(i)*xarm(i));
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
zwt = 0; %Vertical distance bewteen the wing ac and the tail ac
h0 = 0;
a = wCL/AOA; %CL_alpha tail
a1 = 0; %CL_alpha wing
r = lwt/(b/2); 
m = zwt/(b/2);
de_dAOA = 0; %Variation de l'angle epsilon en fonction de l'angle d'attaque
X_np = h0 + V_T*a1/a*(1-(de_dAOA));

%% Aerodynamic center
D_ac = 0.26*(M-0.4)^2.5; %Delta X_ac ; aerodynamic center
X_c4 = 1/4*wMAC; %position of the quarter-chord
X_ac = X_c4 + D_ac*sqrt(S);

%% Static margin
k = X_np - X_cg;
%K = -dC_m/dC_Lw;

% -dC_m/dalpha ; static margin is the difference between the position of 
% the neutral point and the position of the cg or the derivative of the 
% coefficient of the total moment ; -C_malpha/C_Lalpha
% 'Certification authorities specify that k >= 0.05

