function [b,S,CL,CD,D,c_root,c_tip,c_AC,x_AC,y_AC,V_fuel] = wing(Mach,Altitude,Mass)

%% Chosen airfoil: NASA SC(2)-0714  -> Optimal lift coefficient: cl = 0.7
%                                   -> Design/divergence Mach #: MD = 0.725
%                                   -> Maximum thickness ratio: t/c = 0.14
coord   = dlmread('SC(2)-0714.txt');
c_upper = coord(1:103,:);
c_lower = coord(104:end,:);
curve   = [c_upper;c_lower(end:-1:1,:)];

%% Data

% Mach     = 0.7;   % [-]
% Mass     = 5080;  % [kg]
% Altitude = 30000; % [feet]

g  = 9.81;   % [m/s²]
W  = Mass*g; % [N] Take off weight (= Lift at cruise)
L  = W;

[V_inf,rho] = speed(Altitude,Mach);     % 30 000 ft = 9 144 m 
rho = 0.4671+(0.4135-0.4671)*0.144;     % [kg/m³]  Density at 30 000 ft
mu  = (1.458+(1.458-1.493)*0.144)*1e-5; % [kg/m.s] Dynamic viscosity 

AR  = 7;   % Aspect ratio = [5.5-8] for ultralight aircraft
           % Aspect ratio = [7-9] for light aircraft
tap = 0.3; % Tapper ratio (lambda = c_tip/c_root) best value

%% cl and AOA depending on airfoil profil (NASA SC(2)-0714)
% Hypothesis: classical fuselage 
                                           
AOA = 0.75*pi/180; % AOA where the drag is minimum or cl/cd is maximum 
cd  = 0.012;
cl  = 0.7;

%% Computation of wings properties
% Sweep angle
MD    = 0.725;
MObj  = 0.9;
sweep = acos(MD/MObj);
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

CL = a*(AOA-alpha_L0);
% CD = 0.017+CL^2/(0.8*pi*AR);


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
Re   = rho*V_inf*c_root/mu;

% Aerodynamic center
y_AC    = double(2/S*int(c*y,0,b/2));
x_AC    = 0.365*c_AC; % L5 - P26 -> !!! Depend on Mach, Sweep and AR 


%% Computation of fuel volume in the wings

A_tip  = trapz(c_tip*curve(:,1),c_tip*curve(:,2));   % [m²] Area of airfoil
A_root = trapz(c_root*curve(:,1),c_root*curve(:,2)); % [m²] Area of airfoil
V_fuel = (A_root+A_tip)*b/2;                         % [m³] volume of wings

%% Plot of airfoil profil
% 
% figure
% plot(coord(1:103,1),coord(1:103,2),coord(104:end,1),coord(104:end,2))
% plot(c_upper(:,1),c_upper(:,2),c_lower(:,1),c_lower(:,2))
% xlabel('x/c [m]')
% ylabel('z/c [m]')
% axis equal

end