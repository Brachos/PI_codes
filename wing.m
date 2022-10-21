clear all;close all;clc;

%% Data
g     = 9.81;    % [m/s²]
m     = 3500;    % [kg]    Take off mass
W     = m*g;     % [N]     Take off weight (= Lift at cruise)
L     = W;
rho   = 0.4671;  % [kg/m³] Density at 30000 feet
Mach  = 0.7;
Vs    = 303.1;   % [m/s]   Sound speed at 30000 feet
V_inf = Vs*Mach; % [m/s]   Freestream velocity

%%
AR  = 7;   % Aspect ratio = [5.5-8] for ultralight aircraft
           % Aspect ratio = [7-9] for light aircraft
tap = 0.3; % Tapper ratio (lambda = c_tip/c_root) best value


% Computation of cl and cp depending on which NACA
% Hypothesis: classical fuselage 

NACA = [2 4 1 2];
% NACA = [3 3 1 0];

% [cl,cp]  = panel_method(NACA,0,V_inf,1); % Flow can be consider inviscid bc
                                           % V_inf is verry high
cl = 0.25; % Classical cruise lift coefficient according to D. Raymer

% Computation of wings surface

S = 2*L/(rho*V_inf^2*cl);
b = sqrt(AR*S);
c_root = 2*S/((1+tap)*b); % Because for trapez, S = (c_tip+c_root)*b/2
c_tip  = tap*c_root;
