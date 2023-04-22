function [] = main_wing()
%MAIN_WING 

%addpath(genpath('LOADS'));
%addpath(genpath('WING'));

%% Structural loads of the wing
run run_struct_loads

save('wing_loads.mat', 'T_wing', 'M_wing')
%save('fuselage.mat', 'Fuselage')
clear 
clc
load('wing_loads.mat')
%load('fuselage.mat')

%% Material properties of aluminium alloy
SF = 1.5;   %safety factor

sigma_max = 503e6/SF; %tensile Yield Strength [N/m]
tau_max = 331e6/SF;  %shear strength [N/m]
%/!\ valeurs à vérifier :
mu_ref = 27.6e9; %shear modulus [N/m] 
mu = 25.5e9;

%% Input for all the functions computing the structure of the wing
MTOW = 4500; %[kg]
nb_str_root_1 = 4;
nb_str_root_2 = 8;
Mach = 0.7;
Altitude = 3000;
AOA = 2.5;

[cell_root, cell_tip, stringers_root, stringers_tip,span,airf_root] = wing_geom(Mach, Altitude, MTOW,AOA,nb_str_root_1,nb_str_root_2);

%% Boom area computation

coeff_boom1 = 1; %same area for the booms and the stringers

for m = 1 : length(M_wing.X)
[boom_root, boom_tip, stringers_root, area_min(m)] = Boom_area(cell_root,cell_tip,stringers_root, M_wing.X(m), M_wing.Y(m), M_wing.Z(m), sigma_max, coeff_boom1, span);
end
%pq j'avais ajouté span? 

Area = max(area_min);

boom_root.Area(1) = Area*coeff_boom1;
boom_root.Area(2) = Area;
stringers_root.Area = Area;


%% Skin thickness computation

for m = 1 : length(M_wing.X)
[thickness(m)] = skin_size(boom_root,boom_tip, stringers_root, stringers_tip, M_wing.X(m), M_wing.Z(m), T_wing.X(m), T_wing.Z(m),airf_root,mu,mu_ref,tau_max);
end

thickness_max = max(thickness);


%% Display
disp('WING:')
disp(['Stringer area: ', num2str(stringers_root.Area * 1e6), ' [mm2]']);
disp(['Skin thickness: ', num2str(thickness_max* 1e3), ' [mm]']);
disp(['Number of stringers: ', num2str(stringers_root.nb)]);

end

