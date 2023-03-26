%% Code based on Gudmundsson book page 788 

function [CL_max] = takeoff(MTOW,runaway)

% Conversion
    pound = 2.20462262; % kg to lbs
    feet = 3.28; % m to ft
    inche = 39.37; % m to in
    
% Data at cruise
    Mach = 0.7; % Cruise velocity
    Altitude = 30000; % [ft] Altitude at cruise
    AOA = 2.5; % [deg] At root
    
    [~,Sw,~,~,CLw,CD,~,~,~,cw_MAC,~,~,Vw_fuel,~,~] = wing(Mach,Altitude,0.95*MTOW,AOA);
    [D_f,a_el,b_el,l_f,V_f]=fuselage_design(MTOW,Vw_fuel);
    [~,Sh,~,~,~,~,~,CLt,~,~,~,~,~] = v_tail(MTOW,D_f,2*b_el,V_c,cw_MAC
    
    CL_cruise = CLw + Sh/Sw*CLt;
    
% Data at sea level
    rho = 1.225;
    g = 9.81;
    
end 