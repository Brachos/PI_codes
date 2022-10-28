clear all
clc
%%
%design of the fuselage of a BWB


%constraints on the design
M_cruise=0.7;%Mach number at cruise, threshold.
Endurance=10.5;%[h], threshold.
%Vol_payload=;%[m^3], volume necessary to store the payload.
V_cruise=300*M_cruise;%[m/s], cruise speed.
Alt=9144;%[m], flight altitude for cruise (30 000 ft).
M_empty=2250;%[kg]
M_sensors=103.419;%[kg]
M_subsystems=100.698;%[kg]
M_mission=200;%[kg]
M_fuel=1000;%first estimation, half the empty weight of the aircraft.
M_payload=M_sensors+M_mission+M_subsystems+M_fuel;
MTOW=4500;%[kg], maximum take off weight.
g=9.81;%[kg*m/s^2], gravity constant.
rho_cruise=1.123;%[kg/m^3], density of the air at cruise altitude.
Lift=MTOW*g;
c_tip_b=2.25;%[m], is equal to the chord at the root of the wing.
b=2;%[m], span of the body.

%statistical data.
%MTOW and empty weight comparable to a jet trainer (resp. 3500kg and 2500kg).
a=0.333;
C=0.41;
L=a*MTOW^C;%[m], length of the fuselage.
c_root_b=L;%length is in fact the chord at the root of the body. 
T=2;%[m], thickness of the fuselage, based on the volume necessary to store the payload.
tau=T/c_root_b;%thickness to chord ratio.


taper=c_tip_b/c_root_b;


%%
%choose the symetric NACA airfoil to match with tau.
NACA=[0 0 2 1];
%get lift and pressure coefficients
[cl,cp]=panel_method(NACA,0,V_cruise,c_root_b);

S_body=2*Lift/(rho*V_cruise^2*cl);

syms c(y)
c(y)=(c_tip_b-c_root_b)*2/b*y+c_root_b;
integrand=c(y)^2;
MAC=2/S_body*int(integrand,0,b/2);%mean aerodynamic chord.

AR=b^2/S_body;

%%
%the Aspect Ratio of the body is limited by two features:
%-provide a bog enough Cl--> use a surface large enough for the body part.
%-airfoil sufficiently thick to welcolme the payload.
%swept angle of the body 
%taper ratio: ratio bewteen the chord at the root (centerline of the body)
%and the chord at the "tip" (=transition from the body to the wing).
%Several parameters have to be taken into account:
%-chord at the root depends on the thickness-->has to be large enough to
%give enough space to store the payload.
%-tip and root's lengths are chosen to provide enough lift.
% ATTENTION : here the chord at the tip of the body is equal to the
% chord at the root of the wing.

