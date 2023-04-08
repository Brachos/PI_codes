%% Code based on Gudmundsson book page 788 

% function [S_TO] = takeoff(MTOW,CL_max,l_cg)
% Infputs
    % MTOW = Maximum Take-Off Weight
    MTOW = 4231;
    CL_max = 2.5;
    l_cg = 4.51;

% Conversion
    pound = 2.20462262; % kg to lbs
    feet = 3.28; % m to ft
    inche = 39.37; % m to in
    
% Data at cruise
    Mach = 0.7; % Cruise velocity
    Altitude = 30000; % [ft] Altitude at cruise
    AOA = 2.5; % [deg] At root
    
    [bw,Sw,~,~,CLw,CD,~,~,~,cw_MAC,~,~,Vw_fuel,sweep,~] = wing(Mach,Altitude,0.95*MTOW,AOA);
    [~,~,b_el,l_f,~] = fuselage_design(MTOW,Vw_fuel);
    [~,Sh,~,~,~,~,~,CLt,~,~,~,~,~] = v_tail(MTOW,2*b_el,cw_MAC,sweep*180/pi,Sw,l_f,l_cg,bw);
    
    CL_cruise = CLw + Sh/Sw*CLt;
    
% Data at sea level
    rho = 1.225; % [kg/m^3] Air density at sea level 
    mu_air = 1.79e-5; % [Ns/m^2] Air dynamic viscosity at sea level
    g = 9.81;
    mu = 0.04; % Friction coefficient (dry asphalt)
    % mu = 0.05; % "     "     "     "  (wet asphalt)
    AR = 7;    % Wing aspect ratio
    tap = 0.3; % Wing tapper ratio
    h_obst = 30/feet;
    W = MTOW*g;

% Data of engine (check the engine)
    BPR = 4.1; % By-pass ratio
    T = 13300; % [N] Thrust at sea level
    T_mean = 0.75*T*((5 + BPR)/(4 + BPR)); %(18-34 pg 806)
    
% Velocity computation
    Vs = sqrt((2*W)/(rho*CL_max*Sw)); % [m/s] Stall velocity
    V_LOF = 1.1*Vs; % [m/s] Take-off velocity
    V_TR = 1.15*Vs; % [m/s] Transition velocity
    V2 = 1.2*Vs; % [m/s] Climb velocity
    
% Distances computation
    CL_TO = 1.5;%To modify
    CD_TO = 0.017 + CL_TO^2/(0.8*pi*AR);
    D_TO = 0.5*rho*Sw*(V_TR/sqrt(2))^2*CD_TO;
    L_TO = 0.5*rho*Sw*(V_TR/sqrt(2))^2*CL_TO;

    a = g*(T-D_TO-mu*(W-L_TO))/W;
    S_G = V_LOF^2/(2*a);
    S_R = 1*V_LOF; % [m] For small aircraft, t = 1s
    
    gamma = asin((T_mean-D_TO)/W); % [rad]
    R = 0.2156*Vs^2; % [m]
    S_TR = R*sin(gamma); % [m]
    h_TR = R*(1-cos(gamma)); % [m]
    
    if h_TR > h_obst
        S_obst = sqrt(R^2-(R-h_obst)^2);
        S_TO = S_G + S_R + S_obst;
    else
        S_C = (h_obst-h_TR)/tan(gamma);
        S_TO = S_G + S_R + S_TR + S_C;
    end
    
   S=(S_G+S_R)*feet;
% end 