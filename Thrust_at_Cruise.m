%% Description
%
% This is a script to compute an aircraft max thrust at cruise.
%
%% Output
% 
% Display the result.
%
%% Version
% 
% Author: Pavel Azangue (paveleric.azanguedongmo@student.uliege.be)
% Creation date: 13-November-2022
% Matlab version: 2021a
%
%% Revision
%
% V1.0 | 25-February-2023 | Pavel Azangue | Creation
% V1.1
%

%% Program
clear;              % clear workspace
clc;                % clear command Window
close all;          % close all figure

%% 1.) Definitions
%% 1.) - Parameter definition
T_sls = 288.15;                       % Temperature at sea level in [K].
T_alt = 229.733;                      % Temperature at cruis altitude in [K].
M_cruise = 0.75;                      % Mach number at cruise
T_sl = 13.3;                          % Takeoff thrust of the Engine (FJ33-5A) in [KN]
P_sls = 101325;                       % Pressure at sea level [N/m^2]
P_alt = 30800;                        % Pressure at cruise altitude [N/m^2]
rho_alt = 0.46707;                    % Density at cruise altitude [kg/m^3]
rho_sls = 1.225;                      % Density at cruise altitude [kg/m^3]
gamma = 1.4;                          % Air adiabatic parameter.
% G = 0.95;                           % Usable power of a jet engine.
% lambda = 4.1;                       % Bypass ratio.

%% 2.) Computing

% G = input('Enter the gas generator function value: ');

lambda = input('Enter the bypass ratio value : ');

% Computation of temeprature and pressure rates

Theta_0 = (T_alt/T_sls)*(1 + 0.2*M_cruise^2);                     % Temperature rate
delta_0 = (P_alt/P_sls)*(1 + 0.2*M_cruise^2)^3.5;                 % Pressure rate
TR = Theta_0*(1 + 0.2*M_cruise^2);                                % Engine Throtle ratio.

%% Calculation of thrust at cruise for a medium BPR of a turbofan Engine

% A = -0.4327 * (P_alt/P_sls)^2 + 1.3855 * (P_alt/P_sls) + 0.0472;
% X = 0.1377 * (P_alt/P_sls)^2 - 0.4374 * (P_alt/P_sls) + 1.3003;
% Z = 0.9106 * (P_alt/P_sls)^2 - 1.7736 * (P_alt/P_sls) + 1.8697;  
% T_cruise = (((1 - (0.454*(1 + lambda))/sqrt(1 + 0.75*lambda)*G)*M_cruise) + (0.6 + (0.13*lambda)/G)*M_cruise^2 )*T_to;
% T_cruise = (((1 - (0.377*(1 + lambda))/sqrt(1 + 0.82*lambda)*G)*M_cruise) + (0.23 + 0.19*sqrt(lambda))*M_cruise^2 )*T_to;
% T_cruise = ((((A - (0.377*(1 + lambda))/sqrt(1 + 0.82*lambda)*G)* Z *(P_alt/P_sls))*M_cruise) + ((0.23 + 0.19*sqrt(lambda))* X *(P_alt/P_sls))*M_cruise^2 )*T_to;

if Theta_0 < TR
    T_cruise = T_sl*delta_0*(1 - 0.49*sqrt(M_cruise));
elseif Theta_0 > TR
    T_cruise = T_sl*delta_0*(1 - 0.49*sqrt(M_cruise) - 3*(Theta_0 - TR)/(1.5 + M_cruise));
end

%% display of results for analysis

disp(TR);
disp(Theta_0);
disp(delta_0);
disp(T_cruise);
