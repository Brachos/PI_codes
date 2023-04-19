clear
close all

addpath(genpath('LOADS'));
addpath(genpath('FUSELAGE'));

%% Structural loads of the fuselage
run run_struct_loads
save('fuselage_loads.mat', 'T_fuselage', 'M_fuselage')
save('fuselage.mat', 'Fuselage')
clear 
clc
load('fuselage_loads.mat')
load('fuselage.mat')

%% Material properties of aluminium alloy
SF = 1.5;   %safety factor

Stringer.max_load = 503e6/SF; %tensile Yield Strength [N/m]
Stringer.tau_max = 331e6/SF;  %shear strength [N/m]

%% computation of the stringer areas and the skin thickness
[Stringer, Fuselage] = cross_section_geometry(Fuselage, Stringer);

Stringer = stringer(Stringer, T_fuselage, M_fuselage);

thickness = skin(Stringer, T_fuselage, M_fuselage, Fuselage);

%% Display
disp('FUSELAGE:')
disp(['Stringer area: ', num2str(Stringer.B * 1e6), ' [mm2]']);
disp(['Skin thickness: ', num2str(thickness* 1e3), ' [mm]']);
