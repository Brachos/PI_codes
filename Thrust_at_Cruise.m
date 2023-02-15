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
% V1.0 | 13-November-2022 | Pavel Azangue | Creation
% V1.1
%

%% Program
clear;              % clear workspace
clc;                % clear command Window
close all;          % close all figure

%% 1.) Definitions
%% 1.) - Parameter definition

M_cruise = 0.75;                      % Mach number at cruise
T_to = 8.21;                          % Takeoff thrust of the Engine (FJ33-5A) in [KN]
P_sls = 101325;                       % Pressure at sea level [N/m^2]
P_alt = 30800;                        % Pressure at cruise altitude [N/m^2]
rho_alt = 0.467;                      % Density at cruise altitude [kg/m^3]
rho_sls = 1.225;                      % Density at cruise altitude [kg/m^3]
% G = 0.95;                           % Usable power of a jet engine.
% lambda = 4.1;                       % Bypass ratio.

%% 2.) Computing

G = input('Enter the gas generator function value: ');
lambda = input('Enter the bypass ratio value : ');

% A = -0.4327 * (P_alt/P_sls)^2 + 1.3855 * (P_alt/P_sls) + 0.0472;
% X = 0.1377 * (P_alt/P_sls)^2 - 0.4374 * (P_alt/P_sls) + 1.3003;
% Z = 0.9106 * (P_alt/P_sls)^2 - 1.7736 * (P_alt/P_sls) + 1.8697;
    
% T_cruise = (((1 - (0.454*(1 + lambda))/sqrt(1 + 0.75*lambda)*G)*M_cruise) + (0.6 + (0.13*lambda)/G)*M_cruise^2 )*T_to;

 T_cruise = (((1 - (0.377*(1 + lambda))/sqrt(1 + 0.82*lambda)*G)*M_cruise) + (0.23 + 0.19*sqrt(lambda))*M_cruise^2 )*T_to;

% T_cruise = ((((A - (0.377*(1 + lambda))/sqrt(1 + 0.82*lambda)*G)* Z *(P_alt/P_sls))*M_cruise) + ((0.23 + 0.19*sqrt(lambda))* X *(P_alt/P_sls))*M_cruise^2 )*T_to;

disp(T_cruise);
