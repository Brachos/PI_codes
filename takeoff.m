%% Code based on Gudmundsson book page 788 

function [TO] = takeoff(MTOW,WING,V_TAIL,DRAG,WEIGHT,FUSELAGE)

% Conversion
    pound = 2.20462262; % kg to lbs
    feet = 3.28084; % m to ft
    inche = 39.37; % m to in
    
% Wing geometry   
    Sw = WING.Sw;
    cw_root = WING.cw_root;
    cw_tip  = WING.cw_tip;
    AR = WING.A;
    [CL_max_TO, CL_max_L, DCL_max_TO, DCL_max_L, fb_ratio, S_flap] = flaps(WING,FUSELAGE);
    
% Lift and drag coefficient at take-off   
    CL_Cruise = WING.CLw + V_TAIL.Sh_tail*V_TAIL.CL_tail/WING.Sw; % Lift at cruise
    CL_TO = CL_Cruise + DCL_max_TO;
    CD_0 = 0.017;
    e = 1.78*(1-0.045*AR^0.68) - 0.64;
    CD_TO = CD_0 + CL_TO^2/(e*pi*AR);
    
    mu = [0.04 0.05 0.02];
    
for i = 1 : 3    
    
% Data at sea level
    rho = 1.225; % [kg/m^3] Air density at sea level 
    % rho =  (1.112 + 0.524*(1.007-1.112)); % [kg/m^3] Air density at 5000 ft
    % rho =  (0.9093 + 0.048*(0.8194-0.9093)); % [kg/m^3] Air density at 10000 ft
    g = 9.81;
    % mu = 0.04; % Friction coefficient brake off (dry concrete)
    % mu = 0.05; % "        "        "        "   (wet concrete)
    % mu = 0.02; % "        "        "        "   (icy concrete)
    
    h_obst = 30/feet;
    W = MTOW*g;

% Data of engine (check the engine)
    BPR = 4.1; % By-pass ratio
    T = 13580; % [N] Thrust at sea level
    % T = 12170; % [N] Thrust at 5000 ft
    % T = 10740; % [N] Thrust at 10000 ft
    T_mean = 0.75*T*((5 + BPR)/(4 + BPR)); %(18-34 pg 806)
    
% Velocity computation
    Vs = sqrt((2*W)/(rho*CL_max_TO*Sw)); % [m/s] Stall velocity
    V_LOF = 1.1*Vs; % [m/s] Take-off velocity
    V_TR = 1.15*Vs; % [m/s] Transition velocity
    V2 = 1.2*Vs; % [m/s] Climb velocity
    
% Distances computation
    D_TO = 0.5*rho*Sw*(V_TR/sqrt(2))^2*CD_TO;
    L_TO = 0.5*rho*Sw*(V_TR/sqrt(2))^2*CL_TO;

    a = g*(T_mean-D_TO-mu(i)*(W-L_TO))/W;
    S_G(i) = V_LOF^2/(2*a);
    S_R(i) = 1*V_LOF; % [m] For small aircraft, t = 1s
    gamma = asin((T_mean-D_TO)/W); % [rad]
    R = 0.2156*Vs^2; % [m]
    S_TR(i) = R*sin(gamma); % [m]
    h_TR = R*(1-cos(gamma)); % [m]
    S_obst(i) = 0;
    S_C(i) = 0;
    
    if h_TR > h_obst
        S_obst(i) = sqrt(R^2-(R-h_obst)^2);
        S_TO(i) = S_G(i) + S_R(i) + S_obst(i);
    else
        S_C(i) = (h_obst-h_TR)/tan(gamma);
        S_TO(i) = S_G(i) + S_R(i) + S_TR(i) + S_C(i);
    end
    
    S_G(i)  = S_G(i)*feet;  % [ft]
    S_R(i)  = S_R(i)*feet;  % [ft]
    S_TR(i) = S_TR(i)*feet; % [ft]
    S_C(i)  = S_C(i)*feet;  % [ft]
    S_obst(i) = S_obst(i)*feet; % [ft]
    S_TO(i) = S_TO(i)*feet; % [ft]
end
    TO = table(S_G,S_R,S_TR,S_C,S_obst,S_TO);
end 