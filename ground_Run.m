%% Description
% 
% This script compute the ground run of our Aircraft for take-off
% performances and engine validation.
%
%% Output
%
% climb rate and and ground run distance values.
%
%% Version
%
% Author: Pavel Azangue (paveleric.azanguedongmo@student.uliege.be)
% Creation date: 17-December-2022
% Matlab version: 2021a
%
%% Revision
%
% V1.0 | 17-December-2022 | Pavel Azangue | Creation
% V1.1
%

%% Program
clear;              % clear workspace
clc;                % clear command Window
close all;          % close all figure

%% 1.) Definitions
%% 1.) - Parameter definition

% here we assume working with ISA at Sea level

pound = 2.20462262; % kg to lbs
feet = 3.28; % m to ft
inche = 39.37; % m to in

MTOW = 4044;                     % aircraft mass at takeoff in [kg]
Mach = 0.7;
Altitude = 30000; % [ft]
aofa = 0.75; % [°]
[speed,rho] = speed(Altitude,Mach);
rho = 0.48;
V_c = speed;
cg_pos = 4.11;% ? revoir absolument !!!!
l_cg = cg_pos;

[bw,Sw,CLw_alpha,CDw_alpha,CLw,CD,D,cw_root,cw_tip,cw_MAC,xw_AC,yw_AC,Vw_fuel,Lambda_LE,c] = wing(Mach,...
    Altitude,0.95*MTOW,aofa);
[D_f_max,a_el,b_el,l_f,V_f]=fuselage_design(MTOW,Vw_fuel);
[S_tail,S_h,S_v,c_root_tail,c_tip_tail, angle, l, C_L, Lambda_T, b_tail, b_v, b_h, W_tail] = v_tail(MTOW,...
    D_f_max,2*b_el,V_c,cw_MAC,Lambda_LE,Sw,l_f,l_cg,bw);


rho = 1.225;                  % density at sea level in [kg/m^3]
g = 9.80665;                  % gravity value in [m/s^2]
% Sw = 11.76;                     % Reference surface in [m^2]
AR = 7;                       % Aspect ratio
h = 1.5;                     % height of the wing above the ground.
% bw = 9.0700;                      % Wingspan in [m]
CL_cruise = CLw+S_h/Sw*C_L;              % lift coefficient at cruise
mu = 0.03;                    % runway friction coefficient
W = MTOW*g;                      % Aircraft weight at takeoff in [N]
% CL_TO = 1.5;                 % Lift coefficient at take-off
% CL_STR = 0.36;                % Lift coefficient at take-off during transition phase.
% CD_min = 0.025;                % Drag coefficient at take-off
R_T = 0.3;                    % Wing taper ratio. 
% K = 0.1;                      % Wing height in ground roll (it could be less)
beta_L = 0.99;                % High-lift correction coefficient during take-off.
beta_D = 0.1;                 % High-lift correction coefficient during take-off.
n = 1.2;                    % Load factor during take-off
h_obst = 30/feet;               % Obstacle height clearance in [m].
T_SL = 8210;                  % Thrust at sea level in [N].
BPR = 4.1;                    % Engine bypass ratio
K = 0.2;
% n = L/W

%% Computing

% Stall velocity computation

V_stall = sqrt((2*W)/(rho*CL_cruise*Sw));

% computation of the Take-off velocity

V_LOF = 1.1*V_stall;

% computation of the climb velocity

V_TR = 1.15*V_stall;

% Average Thrust during take-off

T_av = 0.75*T_SL*((5 + BPR)/(4 + BPR));

%% Ground effect on Lift and Drag coefficient

% k = (h/b) = 0.1; 

% tapered-wing correction coefficients during take-off.

delta_L = 1 - 2.25*(R_T^0.00273 - 0.997)*(AR^0.717 + 13.6);     

delta_D = 1 - 0.157*(R_T^0.775 - 0.373)*(AR^0.417 - 1.27);    

% High-lift correction coefficient during take-off.

% beta_L = 1 + (0.269 * CL_TO^1.45)/(k^1.12*AR^3.18);

% beta_D = 1 + (0.0361 * CL_TO^1.21)/(k^1.19*AR^1.51);

% Phi_L = (1/beta_L)*(1 + ...
% delta_L*(288*(k^0.787)/AR^0.882)*exp(-9.14*k^0.327));

% Phi_D = beta_D*(1 - delta_D*exp(-4.74*k^0.814) - k^2*exp(-3.88*k^0.758));


% Horizontal segment

q = 0.5*rho*(V_LOF/sqrt(2))^2;              % Dynamic pressure at take-off.

CL_max = 1.5;                               % W/(q*S)*Phi_L;
C_D = 0.1449;                               % (CD_min + k*C_L^2)*Phi_D;

% K_T = (T_SL/W) - mu;
% K_A = ((0.5*rho)/(W/S))*((mu*0.8*CL_max) - C_D - K*(0.8*CL_max)^2);

S_G = 0.5*(V_LOF^2/g)*(1/(6 - (q*Sw*C_D)/W) - mu*(1 - (q*Sw*0.8*CL_max)/W));   % Assumed 70% of Engine Thrust at take-off

% S_G = (1/(2*g*K_A))*log((K_T + K_A*(V_LOF)^2)/K_T);

% Rotation Segment

S_R = 1*V_LOF;                                % For small aircraft the rotation time only take about 1s.

% Transition segment

R = V_TR^2/(0.2*g);                       % Rotation radius during transition phase.

Gama_cl = asind((T_av/W) - (C_D/(0.9*CL_max)));      % Climb angle in  [deg].

S_TR = R*sind(Gama_cl);                            % Transition distance in [m].

h_TR = R*(1 - cosd(Gama_cl));                        % Transition height in [m].

Sc = (h_obst - h_TR)/tand(Gama_cl);                  % Climb distance over an Obstacle.

% S_TR = sqrt(R^2 - (R - h_TR)^2);

% Final Ground Run

S_final = S_G + S_R + S_TR + Sc;

% Display

disp(V_LOF)
disp(h_TR)
disp(S_TR)
disp(S_final)





