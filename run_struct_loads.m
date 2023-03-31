clear

% This script can be run in order to compute the aerodynamic loads as well
% as the structural loads. Most of the parameters are defined in the script
% Input.
%% DATA
g = 9.81; %gravity [m/s²]

run Input
%% Empennage 

Data_empennage=readtable('NACA_visc');  %xfoil computations

AoA_empennage=str2double(Data_empennage.alpha);
CL_empennage=str2double(Data_empennage.CL);
CD_empennage=str2double(Data_empennage.CD);
CM_empennage=str2double(Data_empennage.CM);

a_empennage = 0.065;
%% Wings 

Data_wings=readtable('SC_visc');    %xfoil computations

AoA_wings=str2double(Data_wings.alpha);
CL_wings=str2double(Data_wings.CL);
CD_wings=str2double(Data_wings.CD);
CM_wings=str2double(Data_wings.CM);

beta   = sqrt(1-0.7^2);
L_beta = atan2(tan(sweep),beta);        % Angle to graphically find x_AC

cl_alpha  = (0.76-0.27)*180/(pi*(1+2)); % [1/rad]
alpha_l0  = -(0.83/cl_alpha - 1.5*pi/180);
alpha_01  = -0.23;                      % See L5 - P30
theta_tip = -2*pi/180;                  % [rad] Twist angle
alpha_L0  = alpha_l0 + alpha_01*theta_tip;
k  = beta*cl_alpha/(2*pi);
a_wings  = 2*pi/(2/(beta*Wing.AR)+sqrt((1/(k*cos(L_beta)))^2+(2/(beta*Wing.AR))^2))/beta;

%% Points from the manoeuvre enveloppe
n = [];
V = [];

%% INITIALISATION:
AoA_envelope=zeros(1,length(n));
%Aerodynamic loads
L_W = zeros(1,length(n));
L_E = zeros(1,length(n));
M_fus = zeros(1,length(n));
F_fin = zeros(1,length(n));
D_B = zeros(1,length(n));
M_tail = zeros(1,length(n));
%Strucutural loads: Fuselage
TY =zeros(length(L),length(n));
TZ =zeros(length(L),length(n));
My =zeros(length(L),length(n));
Mz =zeros(length(L),length(n));
Mx =zeros(length(L),length(n));
TX=zeros(length(L),length(n));
%Strucutural loads: wing
Txw=zeros(1,length(n));
Tyw=zeros(1,length(n));
Tzw=zeros(1,length(n));
Mxw=zeros(1,length(n));
Mzw=zeros(1,length(n));
Myw=zeros(1,length(n));

%% COMPUTATION
for i = 1 : length(n)
    Flight.V = V(i);                        %Flight speed [m/s]
    Flight.Mach = Flight.V/speed(30000,1);  %Flight Mach number
    Flight.n = n(i);
    AOA = -25:0.1:25;
    tol = 10^-3;
    error = 1;
    k = 0;
    while error > tol
        k = k + 1;
        Flight.aoa = AOA(k);                %Flight angle of attack [rad]
        L_tot=1/2*Flight.rho*Flight.V^2*Wing.S*a_wings*sind(Flight.aoa+Wing.aoa_fuselage) + 1/2*Flight.rho*Flight.V^2*Empennage.S*a_empennage*sind(Flight.aoa+Empennage.aoa_fuselage) ;

        error = abs(L_tot - n(i)*Aircraft.W);
    end
    AoA_envelope(i) = Flight.aoa;
    if CD_wings(AOA == Flight.aoa)
            Wing.C_L  = a_wings * sind(Flight.aoa + Wing.aoa_fuselage);
            Wing.C_D  = CD_wings(AoA==Flight.aoa) + Wing.C_L^2 / (pi * Wing.AR );
            Wing.C_M  = CM(AoA==Flight.aoa);
            Empennage.C_L = a_empennage * sind(Flight.aoa + Empennage.aoa_fuselage);
            Empennage.C_D = CD(AoA==Flight.aoa) + Empennage.C_L^2 / (pi * Empennage.AR);
            Empennage.C_M =  CM(AoA==Flight.aoa);
    end
    
    %%aerodynamic loads
    [L_W(i), L_E(i), M_fus(i), F_fin(i), D_B(i), M_tail(i)] = aerodynamic_loads(Aircraft, Wing, Empennage, Fin, Flight);
    
    %%Strucutural loads: Fuselage
    [TX(i),TY(i),TZ(i),Mx(i),My(i),Mz(i)] = struct_loads(Fuselage,Wing,Tail,Engine,Sensors,Rear_land_gear,Payload,Empennage,Flight.aoa,n(i), L_E, F_fin, M_fus); 
    
    %%Strucutural loads: Wing 
    Wing_loading.L = L_W(i);
    Wing_loading.D= 1/2 * Flight.rho * Wing.C_D * Flight.V^2 * Wing.S;
    Wing_loading.T = Wing.T;
    AoA_w = (AoA_envelope(i) + Wing.aoa_fuselage) * pi / 180 ;
    y_CG = 0.5;% meter TO BE CHECKED
    y_AC = 0.5;% From P.Raymer TO BE CHECKED
    Mom_wing= 0.5 * Flight.rho * Flight.V^2 * Wing.S * Wing.MAC * Wing.C_M;
    
    [Txw(i),Tyw(i),Tzw(i),Mxw(i),Myw(i),Mzw(i)] = wing_load(Aircraft.W, Wing_loading, AoA_w, Wing.aoa_fuselage * pi / 180, y_cg, y_ac, Mom_wing, n(i), Wing);
end

%% FINAL VALUES
% Fuselage
T_fuselage.X = TX;
T_fuselage.Y = TY;
T_fuselage.Z = TZ;
M_fuselage.X = Mx;
M_fuselage.Y = My;
M_fuselage.Z = Mz;

%Wing root
T_wing.X = Txw;
T_wing.Y = Tyw;
T_wing.Z = Tzw;
M_wing.X = Mxw;
M_wing.Y = Myw;
M_wing.Z = Mzw;