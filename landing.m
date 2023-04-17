%% Code based on Gudmundsson book page 788 

function [S_TO] = landing(MTOW,CL_max,l_cg)
% Inputs
    % MTOW = Maximum Take-Off Weight
    % MTOW = 4243;
    % CL_max = 2.5;
    % l_cg = 4.39;

% Conversion
    pound = 2.20462262; % kg to lbs
    feet = 3.28084; % m to ft
    inche = 39.37; % m to in
    
% Data at cruise
    Mach = 0.7; % Cruise velocity
    Altitude = 30000; % [ft] Altitude at cruise
    AOA = 2.5; % [deg] At root
    
    [bw,Sw,~,~,CLw,CD,~,cw_root,cw_tip,cw_MAC,~,~,Vw_fuel,sweep,~] = wing(Mach,Altitude,0.95*MTOW,AOA);
    [~,~,b_el,l_f,~] = fuselage_design(MTOW,Vw_fuel);
    [~,Sh,~,~,~,~,~,CLt,~,~,~,~,~] = v_tail(MTOW,2*b_el,cw_MAC,sweep*180/pi,Sw,l_f,l_cg,bw);
    
    CL_cruise = CLw + Sh/Sw*CLt;
    
% Data at sea level
    rho = 1.225; % [kg/m^3] Air density at sea level 
    mu_air = 1.79e-5; % [Ns/m^2] Air dynamic viscosity at sea level
    g = 9.81;
    mu = 0.4;  % Friction coefficient brake on (dry concrete)
    % mu = 0.22; % "        "        "        "  (wet concrete)
    % mu = 0.08; % "        "        "        "  (icy concrete)
    
    AR = 7;    % Wing aspect ratio
    tap = 0.3; % Wing tapper ratio
    h_obst = 30/feet;
    W = MTOW*g;

% Data of engine (check the engine)
    BPR = 4.1; % By-pass ratio
    T = 13580; % [N] Thrust at sea level
    T_mean = 0.75*T*((5 + BPR)/(4 + BPR)); %(18-34 pg 806)
    
    