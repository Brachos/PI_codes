%% Code based on Gudmundsson book page 915

function [S_A,S_F,S_FR,S_BR,S_LDG] = landing(MTOW,WING,V_TAIL,DRAG,WEIGHT,FUSELAGE)

% Conversion
    pound = 2.20462262; % kg to lbs
    feet = 3.28084; % m to ft
    inche = 39.37; % m to in
    
% Data at cruise
    Mach = 0.7; % Cruise velocity
    Altitude = 30000; % [ft] Altitude at cruise
    AOA = 2.5; % [deg] At root
    
%     [bw,Sw,~,~,CLw,CD,~,cw_root,cw_tip,cw_MAC,~,~,Vw_fuel,sweep,~] = wing(Mach,Altitude,0.95*MTOW,AOA);
%     [~,~,b_el,l_f,~] = fuselage_design(MTOW,Vw_fuel);
%     [~,Sh,~,~,~,~,~,CLt,~,~,~,~,~] = v_tail(MTOW,2*b_el,cw_MAC,sweep*180/pi,Sw,l_f,l_cg,bw);
    
    Sw = WING.Sw;
    cw_root = WING.cw_root;
    cw_tip  = WING.cw_tip;
    [CL_max_TO, CL_max_L, fb_ratio, S_flap] = flaps(WING,FUSELAGE);
    CL_TO = WING.CLw + V_TAIL.Sh_tail*V_TAIL.CL_tail/WING.Sw; % Lift at cruise
    CD_TO = DRAG.CDTotal;
    
% Data at sea level
    rho = 1.225; % [kg/m^3] Air density at sea level 
    mu_air = 1.79e-5; % [Ns/m^2] Air dynamic viscosity at sea level
    g = 9.81;
    mu = 0.4;   % Friction coefficient brake on (dry concrete)
    % mu = 0.225; % "        "        "        "  (wet concrete)
    % mu = 0.08;  % "        "        "        "  (icy concrete)
    
    AR = 7;    % Wing aspect ratio
    tap = 0.3; % Wing tapper ratio
    h_obst = 50/feet;
    W = WEIGHT.W_empty*g;

% Data of engine (check the engine)
    BPR = 4.1; % By-pass ratio
    T = 13580; % [N] Thrust at sea level
    T_mean = 0.75*T*((5 + BPR)/(4 + BPR)); %(18-34 pg 806)
    
% Velocity computation (using GA aircraft pg 919 Gudmundsson)
    Vs = sqrt((2*W)/(rho*CL_max_L*Sw)); % [m/s] Stall velocity
    V_REF = 1.3*Vs;
    V_FR  = 1.3*Vs;
    V_TD  = 1.1*Vs;
    V_BR  = 1.1*Vs;
    
% Distances computation
    gamma = 3*pi/180; % Approach angle
    h_F = 0.1512*Vs^2*(1-cos(gamma))
    S_A = (h_obst-h_F)/tan(gamma);
    S_F = 0.1512*Vs^2*sin(gamma);
    S_FR = 1*V_TD; % Free-roll distance (1 sec for small aircraft)
    
    T_LDG = -0.4*T;
    CD_0 = 0.017;
    e = 1.78*(1-0.045*AR^0.68) - 0.64;
    CL_LDG = WING.CLw + V_TAIL.Sh_tail*V_TAIL.CL_tail/WING.Sw; % Lift at cruise
    CD_LDG = DRAG.CDTotal;
    L_LDG = 0.5*rho*CL_LDG*Sw*(V_BR/sqrt(2))^2;
    D_LDG = 0.5*rho*CD_LDG*Sw*(V_BR/sqrt(2))^2;
    
    a = g/W*(T_LDG-D_LDG-mu*(W-L_LDG))
    S_BR = -V_BR^2/a;
    
    S_LDG = S_A + S_F + S_FR + S_BR;    
    
    
%     CL_TO = 0.8;%To modify
%     CD_TO = 0.017 + CL_TO^2/(0.8*pi*AR);
%     D_TO = 0.5*rho*Sw*(V_TR/sqrt(2))^2*CD_TO;
%     L_TO = 0.5*rho*Sw*(V_TR/sqrt(2))^2*CL_TO;
    
end  