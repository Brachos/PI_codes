%% Parameters
b = 10.25; %wingspan [m]
c_ = 1.463; %standard mean chord (smc) [m]
c__ = 0; %mean aerodynamic chord (mac)
S_T = 0; %Tailplane area
S_F = 0; %Fin area
l_T = 0; %Tail moment arm
l_F = 0; %Fin moment arm
le = 0; %length of the engine
L = 0; %length of the fuselage
Wt = 0; %minimum weight (or minimum mass)
Mt = zeros(3,1); %vector of the total moment 3 different directions
M = 0.7; %Mach number
x_wLE = 0; %position of the leading edge of the wMAC
x_vLE = 0; %position of the leading edge of the vMAC
x_e = 0; %position of the engine inlet
Nelem = 0; % number of differents elements, of different mass
           % (1.Wing;2.Fuselage;3.Tail;4.Engines;5.Landing gears;
           % 6.Payload;7.Fuel?)
wrc = 2.25; %wing root chord [m]
wtc = 0.68; %wing tip chord [m]
wTR = wtc/wrc; %wing taper ratio
vrc = 0; %v-tail root chord
vtc = 0; %v-tail tip chord 
vTR = 0; %v-tail taper ratio
%% Wing Main Aerodynamic Chord (MAC)
wMAC = wrc*(2/3)*((1+wTR+wTR^2)/(1+wTR));
%% V-Tail Main Aerodynamic Chord
vMac = vrc*(2/3)*((1+vTR+vTR^2)/(1+vTR));

%% 
S = b*c_; %wing area S=15 [m]
AR = b^2/S; %Aspect Ratio (from 1.5 to 18)
V_T = S_T*l_T/S*c__; %Tail volume ratio
V_F = S_F*l_F/S*c__; %fin volume ratio

%% Weight
W = []; %vector of all the different weights (or mass)
MTOW = sum(W); %Maximum Take-Off Weight

%% Center of gravity
xcg_w = 0.4*wMAC;  %for the wing [35%wMAC;42%wMAC] [m]
xcg_e = 0.37*le; %for the engine [30%le;45%le] [m]
xcg_f = 0.44*L; %for the fuselage [40%L;48%L] [m]
xcg_v = 0.3*vMAC; %for the v-tail [m]
xcg_l= 0; %for the landing gears
xcg_p = 0; %for the payload
thick = 0; %thickness of the plane
zcg_w = 0.07*thick; %for the wing [0.05*thick;0.10*thick] [m]
xarm = [xcg_w+x_wLE;xf;xcg_v+x_vLE;wcg_e+x_e;xcg_l;xcg_p;0];
yarm = [0;0;0;0;0;0;0]; %symetric
zarm = [zcg_w;0;0;0;0;0;0];
% vector of all the different arms corresponding to the different
% weights ; arm = horizontal dimension from the nose of the plane
% to the cg of the mass ; origin : nose
cgT = zeros(3,1); %coordinates of the center of gravity maximum weight
cgt = zeros(3,1); %coordinates of the center of gravity minimum weight
X_cg=0;
for i=1:Nelem
    for j=1:3
        Mt(j) = Mt(j) +(W(i)*xarm(i));
    end
end
for i=1:3
    cgT(i) = Mt(i)/MTOW; 
    cgt(i) = Mt(i)/Wt;
    % Total moment/Total weight ; we can compute a range for the cg with
    % the lower weight and the greater weight different coordinates of the 
    % cg for the max weight
end

%% Neutral point
X_np = 0;

%% Aerodynamic center
D_ac = 0.26*(M-0.4)^2.5; %Delta X_ac ; aerodynamic center
X_c4 = 0; %position of the quarter-chord
X_ac = X_c4 + D_ac*sqrt(S);

%% Static margin
k = X_np - X_cg; 
% -dC_m/dalpha ; static margin is the difference between the position of 
% the neutral point and the position of the cg or the derivative of the 
% coefficient of the total moment ; -C_malpha/C_Lalpha
% 'Certification authorities specify that k >= 0.05

