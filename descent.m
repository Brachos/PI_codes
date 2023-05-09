%% Code based on Gudmundsson book page 915

function [DESCENT] = descent(MTOW,WING,WEIGHT)

% Conversion
    pound = 2.20462262; % kg to lbs
    feet = 3.28084; % m to ft
    inche = 39.37; % m to in
% Data at cruise
    Mach = 0.7; % Cruise velocity
    Altitude = 30000; % [ft] Altitude at cruise
    
    Sw = WING.Sw;

% Data at 30 000 ft
    h = 30000/feet;
    rho = 0.48; % [kg/m^3] Air density at 30 000 ft 
    g = 9.81;
    
    AR = 7;    % Wing aspect ratio
    tap = 0.3; % Wing tapper ratio

    W_TO = MTOW*g;
    W_LDG = WEIGHT.W_empty*g;
    Sw = WING.Sw;
    CD_0 = 0.019;
    e = 1.78*(1-0.045*AR^0.68) - 0.64;
    
    % Case of engine failure
    
    V_BR_TO = sqrt(2*W_TO/(rho*Sw)*sqrt(1/(CD_0*pi*e*AR))); % [m/s]
    V_BR_LDG = sqrt(2*W_LDG/(rho*Sw)*sqrt(1/(CD_0*pi*e*AR))); % [m/s]
    L_D_max = 1/sqrt(4*CD_0/(e*pi*AR));
    R_glide = h*L_D_max; % [m]
    ROD_TO = V_BR_TO*sin(atan2(h,R_glide)); % [m/s]
    ROD_LDG = V_BR_LDG*sin(atan2(h,R_glide)); % [m/s]
    
    V_BR_TO = feet*V_BR_TO;
    V_BR_LDG = feet*V_BR_LDG;
    R_glide = feet*R_glide;
    ROD_TO = feet*ROD_TO;
    ROD_LDG = feet*ROD_LDG;
    
    DESCENT = table(V_BR_TO,V_BR_LDG,L_D_max,R_glide,ROD_TO,ROD_LDG);
    
    % Normal case
    
    
end
    