function [b,S,CL_alpha,CD_alpha,CL,CD,D,c_root,c_tip,c_AC,x_AC,y_AC,V_fuel,sweep,c,alpha_L0, tap] = wing(Mach,Altitude,Mass,AOA)

%% Chosen airfoil: NASA SC(2)-0714  -> Optimal lift coefficient: cl = 0.7
%                                   -> Design/divergence Mach #: MD = 0.725
%                                   -> Maximum thickness ratio: t/c = 0.14
coord   = dlmread('SC(2)-0714.txt');
c_upper = coord(1:103,:);
c_lower = coord(104:end,:);
curve   = [c_upper;c_lower(end:-1:1,:)];

%% Data

% AOA      = 2.5;  % [deg]
% Mach     = 0.7;   % [-]
% Mass     = 4471;  % [kg]
% Altitude = 30000; % [feet]

g  = 9.81;   % [m/s^2]
L  = Mass*g; % [N] Take off weight (= Lift at cruise)

[V_inf,rho] = speed(Altitude,Mach);     % 30 000 ft = 9 144 m 
rho = 0.4671+(0.4135-0.4671)*0.144;     % [kg/m^3] Density at 30 000 ft
mu  = (1.458+(1.458-1.493)*0.144)*1e-5; % [kg/m.s] Dynamic viscosity 

AR  = 7;   % Aspect ratio = [5.5-8] for ultralight aircraft
           % Aspect ratio = [7-9] for light aircraft
tap = 0.3; % Tapper ratio (lambda = c_tip/c_root) best value
                                           
AOA = AOA*pi/180; % AOA where the drag is minimum or cl/cd is maximum 

%% Computation of wings properties
% Sweep angle
% MD    = 0.725;
% MObj  = 0.8;
% sweep = acos(MD/MObj);
sweep  = 15*pi/180; % [rad] Sweep angle 
beta   = sqrt(1-Mach^2);
L_beta = atan(tan(sweep)/beta); % Angle to graphically find x_AC

% CL and CD
cl_alpha  = (0.76-0.27)*180/(pi*(1+2)); % [1/rad]
alpha_l0  = -(0.83/cl_alpha - 1.5*pi/180);
alpha_01  = -0.23;     % See L5 - P30
theta_tip = -2*pi/180; % [rad] Twist angle
alpha_L0  = alpha_l0 + alpha_01*theta_tip;
k  = beta*cl_alpha/(2*pi);
a  = 2*pi/(2/(beta*AR)+sqrt((1/(k*cos(L_beta)))^2+(2/(beta*AR))^2))/beta;
CL_alpha = a;

CL = a*(AOA-alpha_L0);
CD = 0.017+CL^2/(0.8*pi*AR); % CD_0 = 0.017 and e = 0.8 for high subsonic 

CD_alpha = 2*a^2*(AOA-alpha_L0)/(0.8*pi*AR);

% Wing surface and span
S  = 2*L/(rho*V_inf^2*CL);
b  = sqrt(AR*S);

% Thrust
D  = 0.5*CD*rho*V_inf^2*S;

% Chord
c_root = 2*S/((1+tap)*b); % Because for trapez, S = (c_tip+c_root)*b/2
c_tip  = tap*c_root;
syms y
c    = (1-2*y/b)*c_root + (2*y/b)*c_tip;
c_AC = double(2/S*int(c^2,0,b/2));

% Aerodynamic center
y_AC    = double(2/S*int(c*y,0,b/2));
x_AC    = 0.365*c_AC; % L5 - P26 -> !!! Depend on Mach, Sweep and AR 

% Flaps design

CL_max_base = 0.95*cos(sweep)*1.75*1.1 % Maximum CL without device
CL_max = 2.5;
DCL_max = CL_max - CL_max_base
S_flapped = S*DCL_max/(1.3*0.9*cosd(10)) % 1.3 because slotted flap

%% Computation of fuel volume in the wings

A_tip  = trapz(c_tip*curve(:,1),c_tip*curve(:,2));   % [m^2] Area of airfoil
A_root = trapz(c_root*curve(:,1),c_root*curve(:,2)); % [m^2] Area of airfoil
V_fuel = (A_root+A_tip)*b/2*0.8^2;                   % [m^2] volume of wings

%% Plot of airfoil profil
 
% figure
% plot(curve(:,1),curve(:,2),'k','LineWidth',1)
% axis equal
% axis off

end