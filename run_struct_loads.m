clear

% This script can be run in order to compute the aerodynamic loads as well
% as the structural loads. Most of the parameters are defined in the script
% Input.
%% DATA
g = 9.81; %gravity [m/s²]
AOA = -50:0.1:50;                       %Angle of attack of the wing [deg]
tol = 1e-1;
    
run Input
%% Empennage 

Data_empennage=readtable('NACA_visc');  %xfoil computations

AoA_empennage=Data_empennage.alpha;
CL_empennage=Data_empennage.CL;
CD_empennage=Data_empennage.CD;
CM_empennage=Data_empennage.CM;

a_empennage = 5.7296;      %from the tail design
%% Wings 

Data_wings=readtable('SC_visc.txt');    %xfoil computations

AoA_wings=Data_wings.alpha;
CL_wings=Data_wings.CL;
CD_wings=Data_wings.CD;
CM_wings=Data_wings.CM;

a_wings  = 6.015050996705918; %from the wing design

%% Points from the manoeuvre enveloppe
n = [4.5 -1.5 3 -1.5 3]; 
V = [198.48 198.48 228.25 228.25 113]; %[m/s]

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
TY =zeros(1,length(n));
TZ =zeros(1,length(n));
My =zeros(1,length(n));
Mz =zeros(1,length(n));
Mx =zeros(1,length(n));
TX =zeros(1,length(n));

%Strucutural loads: wing
Txw=zeros(1,length(n));
Tyw=zeros(1,length(n));
Tzw=zeros(1,length(n));
Mxw=zeros(1,length(n));
Mzw=zeros(1,length(n));
Myw=zeros(1,length(n));
%% Angle of attack
 for i = 1: length(n)                        %Flight speed [m/s]    
     error = 1;
     k = 0;
     while error > tol
         k = k + 1;
         L_tot=1/2*Flight.rho*V(i)^2*Wing.S*a_wings*sind(AOA(k)) + 1/2*Flight.rho*V(i)^2*Empennage.S*a_empennage*sind(AOA(k) - Wing.aoa_fuselage*180/pi + Empennage.aoa_fuselage*180/pi) ;
         
         error = abs((L_tot - n(i)*Aircraft.W)/(n(i)*Aircraft.W));
     end
     
     AoA_envelope(i) = AOA(k);
 end
%% COMPUTATION
for i = 1 : length(n)
    Flight.V = V(i);                        %Flight speed [m/s]
    Flight.Mach = Flight.V/speed(30000,1);  %Flight Mach number
    Flight.n = n(i);
    
    Flight.aoa = AoA_envelope(i);
    if abs(Flight.aoa) < 10
         Wing.C_L  = a_wings * sind(Flight.aoa);
         Wing.C_D  = interp1(AoA_wings, CD_wings, Flight.aoa) + Wing.C_L^2 / (pi * Wing.AR );
         Wing.C_M  = interp1(AoA_wings, CM_wings, Flight.aoa);
    else
        Wing.C_L = a_wings * sind(Flight.aoa);
        Wing.C_D = (CD_wings(AoA_wings == 10) - CD_wings(AoA_wings == 0))/10*Flight.aoa + Wing.C_L^2 / (pi * Wing.AR );
        Wing.C_M = (CM_wings(AoA_wings == 10) - CM_wings(AoA_wings == 0))/10*Flight.aoa;
    
    end
    
    if abs(Flight.aoa -Wing.aoa_fuselage*180/pi + Empennage.aoa_fuselage*180/pi) < 5
        Empennage.C_L = a_empennage * sind(Flight.aoa -Wing.aoa_fuselage*180/pi + Empennage.aoa_fuselage*180/pi);
        Empennage.C_D = interp1(AoA_empennage, CD_empennage, Flight.aoa-Wing.aoa_fuselage*180/pi + Empennage.aoa_fuselage*180/pi) + Empennage.C_L^2 / (pi * Empennage.AR);
        Empennage.C_M = interp1(AoA_empennage, CM_empennage, Flight.aoa-Wing.aoa_fuselage*180/pi + Empennage.aoa_fuselage*180/pi);
    else
        if Flight.aoa -Wing.aoa_fuselage*180/pi + Empennage.aoa_fuselage*180/pi > 5
            AOA = Flight.aoa-Wing.aoa_fuselage*180/pi + Empennage.aoa_fuselage*180/pi;
            Empennage.C_L = a_empennage * sind(Flight.aoa -Wing.aoa_fuselage*180/pi + Empennage.aoa_fuselage*180/pi);
            Empennage.C_D = (CD_empennage(AoA_empennage == 5) - CD_empennage(AoA_empennage == 0))/5*AOA + Empennage.C_L^2 / (pi * Empennage.AR);
            Empennage.C_M = (CM_empennage(AoA_empennage == 5) - CM_empennage(AoA_empennage == 0))/5*AOA;
        else
            AOA = Flight.aoa-Wing.aoa_fuselage*180/pi + Empennage.aoa_fuselage*180/pi;
            Empennage.C_L = a_empennage * sind(Flight.aoa -Wing.aoa_fuselage*180/pi + Empennage.aoa_fuselage*180/pi);
            Empennage.C_D = (CD_empennage(AoA_empennage == 0) - CD_empennage(AoA_empennage == -5))/5*(-AOA) + Empennage.C_L^2 / (pi * Empennage.AR);
            Empennage.C_M = (CM_empennage(AoA_empennage == 0) - CM_empennage(AoA_empennage == -5))/5*(-AOA);
        end
    end
    
    Flight.aoa = deg2rad(Flight.aoa);
    %%aerodynamic loads
    [L_W(i), L_E(i), M_fus(i), F_fin(i), D_B(i), M_tail(i)] = aerodynamic_loads(Aircraft, Wing, Empennage, Fin, Flight);
    
    %%Strucutural loads: Fuselage
    [TX(i), TY(i), TZ(i), Mx(i), My(i), Mz(i)] = struct_loads(Fuselage,Wing,Tail,Engine,Sensors,Rear_land_gear,Payload,Empennage,Fin,Flight.aoa,Flight.n, L_E(i), F_fin(i), M_fus(i)); 
    
    %%Strucutural loads: Wing 
    Wing_loading.L = L_W(i);
    Wing_loading.D = 1/2 * Flight.rho * Wing.C_D * Flight.V^2 * Wing.S;
    y_ac = 1.7168;
    y_CG = 1.2650;% meter TO BE CHECKED
    Mom_wing= 0.5 * Flight.rho * Flight.V^2 * Wing.S * Wing.MAC * Wing.C_M;
    
    [Txw(i),Tyw(i),Tzw(i),Mxw(i),Myw(i),Mzw(i)] = wing_load(Wing.W/2, Wing_loading, Flight.aoa, y_CG, y_ac, Mom_wing, Flight.n, Wing);
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