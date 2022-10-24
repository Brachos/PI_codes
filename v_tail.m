%% Données
AR= 7;
SW = 10.76; %m²
bW = 10.8; %m
c_root = 1.9; %m
c_tip = 0.57; %m
cVT_init = 0.06;  %between 0.02-0.09
cHT_init = 0.8; %between 0.5-1
LVT = 1; %vertical distance [m]
LHT = 8; %horizontal distance [m]
%coucou
%yolo
%% Calculations
fun = @(x) (c_root + (c_tip-c_root)/bW*2*x).^2;
CW_bar = 2/SW * integral(fun,0,bW/2);
SVT_init = cVT_init*bW*SW/LVT;
SHT_init = cHT_init*CW_bar*SW/LHT;
dihedral_angle_init = atan(sqrt(SVT_init/SHT_init));
fprintf("L'angle initial dihedre est de %d degrés\n", dihedral_angle_init*180/2/pi)
% choix a priori d'un profil NACA 0012
