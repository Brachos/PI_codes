clear all
clc
%% 
%conditions of flight
M_cruise=0.7;%Mach number at cruise, threshold.
Endurance=10.5;%[h], threshold.
Alt=9144;%[m], flight altitude for cruise (30 000 ft).
V_cruise=300*M_cruise;%[m/s], cruise speed.
g=9.81;%[kg*m/s^2], gravity constant.
rho_cruise=1.123;%[kg/m^3], density of the air at cruise altitude.
%%
%mass of the payload
M_sensors=103.419;%[kg]
M_subsystems=100.698;%[kg]
M_mission=200;%[kg], water + radio
M_fuel=1000;%[kg]
M_payload=M_sensors+M_mission+M_subsystems+M_fuel;
M_empty=2250;%[kg]
%MTOW=M_payload+M_empty;
MTOW=4471;%[kg], first estimation.
%%
%volume of the payload
% V_sensors=;
% V_subsystems=;
V_mission=0.15;%[m^3]
f=0;%fraction of fuel in the fuselage
% V_fuel=;
% V_landing_gear=;
%V_int_fuselage=V_sensors+V_subsystems+V_mission+V_fuel+V_landing_gear;
%%
%first estimation using statistical relations.
%MTOW and empty weight comparable to a jet trainer (resp. 3500kg and 2500kg).
a=0.333;
C=0.41;
L_first_guess=a*MTOW^C;%[m], length of the fuselage.
%%
%Volume of the fuselage, kept constant through the iterations. 
%IMPORTANT: fuselage is first considered to be cylindrical and will later
%be streamlined.
%The volume used here corresponds to the volume necessary to store all the
%payload+the subsystems+the sensors+the landing gears+fuel.
V_f=3;%[m^3]
f=7;%optimum fineness ratio of the fuselage to minimize drag (for a constant volume):[6;8].
%first estimation of the diameter
D_m_first_guess=2*sqrt(V_f/(pi*L_first_guess));
%second estimation
L_sec_guess=f*D_m_first_guess;
D_sec_guess=2*sqrt(V_f/(pi*L_sec_guess));
%iterations to get the final length and diameter of the fuselage.
L_f=zeros(2,1);
D_f=zeros(2,1);
L_f(1)=L_first_guess;
L_f(2)=L_sec_guess;
D_f(1)=D_m_first_guess;
D_f(2)=D_sec_guess;
i=2;
while abs(L_f(i)-L_f(i-1))>0.01 && abs(D_f(i)-D_f(i-1))>0.001
    L_f(i+1)=f*D_f(i);
    D_f(i+1)=2*sqrt(V_f/(pi*L_f(i+1)));
    i=i+1;
end
fprintf('The final length of the fuselage is equal to %d m \n',L_f(i));
fprintf('The final equivalent diameter of the fuselage is equal to %d m \n',D_f(i));
a=D_f(i)/sqrt(2);
b=a/2;
A_f=a*b*pi;
V_f=A_f*L_f(i);
fprintf('Dimensions of the elliptical cross-section: a=%d m and b=%d m \n',a,b);
fprintf('The area of the elliptical cross-section is equal to %d m^2 \n',A_f);
fprintf('Volume of the fuselage: V=%d m^3 \n',V_f);