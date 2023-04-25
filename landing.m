%% Code based on Gudmundsson book page 915

function [LDG] = landing(MTOW,WING,V_TAIL,DRAG,WEIGHT,FUSELAGE)

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
    
% Lift and drag coefficient at landing
    CL_Cruise = WING.CLw + V_TAIL.Sh_tail*V_TAIL.CL_tail/WING.Sw; % Lift at cruise
    CL_LDG = CL_Cruise + DCL_max_L;
    CD_0 = 0.017;
    e = 1.78*(1-0.045*AR^0.68) - 0.64;
    CD_LDG = CD_0 + CL_LDG^2/(e*pi*AR);

    mu = [0.4 0.225 0.08];
    
for i = 1 : 3
    
% Data at sea level
    rho = 1.225; % [kg/m^3] Air density at sea level     
    % rho =  (1.112 + 0.524*(1.007-1.112)); % [kg/m^3] Air density at 5000 ft
    % rho =  (0.9093 + 0.048*(0.8194-0.9093)); % [kg/m^3] Air density at 10000 ft
    % rho =  (0.8194 + 0.572*(0.7364-0.8194)); % [kg/m^3] Air density at 15000 ft
    % rho =  (0.6601 + 0.096*(0.5900-0.6601)); % [kg/m^3] Air density at 20000 ft
    
    g = 9.81;
    % mu = 0.4;   % Friction coefficient brake on (dry concrete)
    % mu = 0.225; % "        "        "        "  (wet concrete)
    % mu = 0.08;  % "        "        "        "  (icy concrete)
    
    h_obst = 50/feet;
    W = WEIGHT.W_empty*g;

% Data of engine (check the engine)
    BPR = 4.1; % By-pass ratio
    T = 13580; % [N] Thrust at sea level
    % T = 12170; % [N] Thrust at 5000 ft
    % T = 10740; % [N] Thrust at 10000 ft
    % T = 9327; % [N] Thrust at 15000 ft
    % T = 7985; % [N] Thrust at 20000 ft
    
% Velocity computation (using GA aircraft pg 919 Gudmundsson)
    Vs = sqrt((2*W)/(rho*CL_max_L*Sw)) % [m/s] Stall velocity
    V_REF = 1.3*Vs;
    V_FR  = 1.3*Vs;
    V_TD  = 1.1*Vs;
    V_BR  = 1.1*Vs;
    
% Distances computation
    gamma = 3*pi/180; % Approach angle
    h_F = 0.1512*Vs^2*(1-cos(gamma))
    S_A(i) = (h_obst-h_F)/tan(gamma);
    S_F(i) = 0.1512*Vs^2*sin(gamma);
    S_FR(i) = 1*V_TD; % Free-roll distance (1 sec for small aircraft)
    
    T_LDG = -0.4*T;
    L_LDG = 0.5*rho*CL_LDG*Sw*(V_BR/sqrt(2))^2;
    D_LDG = 0.5*rho*CD_LDG*Sw*(V_BR/sqrt(2))^2;
    
    a = g/W*(T_LDG-D_LDG-mu(i)*(W-L_LDG));
    S_BR(i) = -V_BR^2/a;
    S_LDG(i) = S_A(i) + S_F(i) + S_FR(i) + S_BR(i);
    
    S_A(i) = S_A(i)*feet; % [ft]
    S_F(i) = S_F(i)*feet; % [ft]
    S_FR(i) = S_FR(i)*feet; % [ft]
    S_BR(i) = S_BR(i)*feet; % [ft]
    S_LDG(i) = S_LDG(i)*feet; % [ft]
end
    LDG = table(S_A,S_F,S_FR,S_BR,S_LDG);
end  