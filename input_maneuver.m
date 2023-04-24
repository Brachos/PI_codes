%%input for the maneuver enveloppe

in.MTOW = 3.3330e+03;%[kg]

%extreme load factor:
n.max = 3;%min(2.1 + 24000/(10000 + in.MTOW),3.8); % KPP19 
n.min = -1.5;% KPP18
in.S = 7.3885; %[m^2]

in.cw_MAC = 1.1267; %mean aero chord [m]
in.C_L_max = 1.6655;   %in cruise
in.C_L_max0 = -2.2; %at landing
in.C_L_alpha_plane = 5.99982543977986; %[1/rad]

Ue.Vb = (20.86-44)/(60000-15000)*15000 + 44; %[ft/s] 36
Ue.Vc = Ue.Vb; %[ft/s]
Ue.Vd = (10.43-22)/(60000-15000)*15000 + 22; %[ft/s]  18

%Ue.Vb = 56.67; %[ft/s]
%Ue.Vc = 41.67; %[ft/s]
%Ue.Vd = 20.83; %[ft/s]

%these values are obtained interpolating the FAR rule page 9 "discrete gust
%static loads"