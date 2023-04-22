function [S_TO,S_LDG] = TO_LDG(WING,V_TAIL,FUSELAGE,MTOW,WEIGHT,DRAG)

% Conversion
    pound = 2.20462262; % kg to lbs
    feet = 3.28084; % m to ft
    inche = 39.37; % m to in
% Data at cruise
    Mach = 0.7; % Cruise velocity
    Altitude = 30000; % [ft] Altitude at cruise
    
    Sw = WING.Sw;
    cw_root = WING.cw_root;
    cw_tip  = WING.cw_tip;
    [CL_max_TO,CL_max_LDG,~,~] = flaps(WING,FUSELAGE);
    CL_TO = WING.CLw + V_TAIL.Sh_tail*V_TAIL.CL_tail/WING.Sw; % Lift at cruise
    CD_TO = DRAG.CDTotal;
    
% Data at sea level
    rho = 1.225; % [kg/m^3] Air density at sea level 
    mu_air = 1.79e-5; % [Ns/m^2] Air dynamic viscosity at sea level
    g = 9.81;
    mu_TO = 0.04; % Friction coefficient brake off (dry concrete)
    % mu_TO = 0.05; % "        "        "        "   (wet concrete)
    % mu_TO = 0.02; % "        "        "        "   (icy concrete)
    mu_LDG = 0.4;   % Friction coefficient brake on (dry concrete)
    % mu_LDG = 0.225; % "        "        "        "  (wet concrete)
    % mu_LDG = 0.08;  % "        "        "        "  (icy concrete)
    
    mu_TO = mu_TO + 0.72*0.017/CL_max_TO;
    
    AR = 7;    % Wing aspect ratio
    tap = 0.3; % Wing tapper ratio
    h_obst = 30/feet;
    W_TO = MTOW*g;
    W_LDG = WEIGHT.W_empty*g;

    BPR = 4.1; % By-pass ratio
    T = 13580; % [N] Thrust at sea level
    T_mean = 0.75*T*((5 + BPR)/(4 + BPR)); %(18-34 pg 806)

    
%% Takeoff    
% Velocity computation
    Vs = sqrt((2*W_TO)/(rho*CL_max_TO*Sw)); % [m/s] Stall velocity
    V_LOF = 1.1*Vs; % [m/s] Take-off velocity
    V2 = 1.2*Vs; % [m/s] Climb velocity
    
% Distance computation
    S_G = (V_LOF^2/(2*g))/(T_mean/W_TO - mu_TO);
    gamma_LOF = 0.9*T_mean/W_TO - 0.3/sqrt(AR);
    S_A = V_LOF^2/(g*sqrt(2)) + h_obst/gamma_LOF;
    S_tot = S_G + S_A;
    S_TO = table(S_tot,S_G,S_A);
    
%% Landing
% Velocity computation
    Vs = sqrt((2*W_LDG)/(rho*CL_max_LDG*Sw)) % [m/s] Stall velocity
    V_REF = 1.2*Vs;
    V_TD = 1.1*Vs;
    
% Dictance computation
    gamma = 3*pi/180; % Approach angle
    h_F = 0.1512*Vs^2*(1-cos(gamma));
    x3 = (h_obst-h_F)/tan(gamma);
    t3 = x3/(V_REF*cos(gamma));
    n = 1.296; % Load factor
    R = V_REF^2/(g*(n-1));
    x2 = R*sin(gamma);
    t2 = gamma*V_REF/(g*(n-1));
    a = 0.6; % For light aircraft with only brake
    xg = V_TD^2/a;
    x_tot = x3 + x2 + xg;
    S_LDG = table(x_tot,x3,x2,xg);
end