function [HE,HT,Hmfg,N_eng,t_ac,CPI,Ceng,Cdev,CFT,Ctool,CMFG,Cqc,Cmat,Ccert,Cpp,cost_per_aircraft] = general_aviation(We,V_max,N,T)
    %% Computation of the number of hours page 49 du livre en version pdf 
    Fcert = 1;
    Fcf = 1;
    Fcomp = 1;
    Fpress = 1; %unpressurized aircraft
    Ftaper = 1; %tappered wing
    HE = 0.0396 * We^(0.791) * V_max^(1.526) * N^(0.183)* Fcomp*Fcert*Fcomp*Fpress; %Engineering hours
    HT = 1.0032 * We^(0.764) * V_max^(0.899) * N^(0.178) * (N/60)^(0.066)*Ftaper*Fcf*Fcomp*Fpress; % Tooling hours
    Hmfg = 9.6613*We^(0.74)*V_max^(0.543)*N^(0.524)*Fcert*Fcf*Fcomp; %Manufacturing hours
    
    %% Number of engineer to developp the aircraft in 5 years and average time to manufacture
    N_eng = HE/(5*48*40) % Number of engineer to developp the aircraft in 5 years
    t_ac = Hmfg/N;
    
    %% Fixed cost and variable cost
    Np = 1; %Number of prototype
    Npp =1; %Number of engine
    year = 2012:1:2022;
    CPI = [1.021 1.015 1.016 1.001 1.013 1.021 1.024 1.018 1.012 1.047 1.086];
    Reng = 92; %Rate engineering labour [dollar/hour]
    Rmfg = 53; %Rate of manufacturing labour [dollar/hour]
    Rtool = 61; %Rate of tooling labour [dollar/hour]
    Reng_evolution = []
    Reng_evolution(1) = 92
    for i=2:length(CPI)
        Reng_evolution(i) = Reng_evolution(i-1)*CPI(i)
    end
    CPI = Reng_evolution(end)/Reng; %CPI 
    Ceng = 2.0969*HE*CPI *Reng; %Total Cost of engineering
    Cdev = 0.06458*We^(0.873)*V_max^(1.89)*Np^(0.326)*CPI*Fcert*Fcf*Fcomp*Fpress; %Total cost of dev
    CFT = 0.009646*We^(1.16)*V_max^(1.3718)*Np^(1.281)*CPI*Fcert; %Total flight test operation
    Ctool = 2.0969*HT*Rtool*CPI;
    CMFG = 2.0969*Hmfg*Rmfg*CPI; % Total cost of manufacturing
    Cqc = 0.13*CMFG*Fcert*Fcomp; % Total cost of quality control
    Cmat = 24.896*We^(0.689)*V_max^(0.624)*N^(0.792)*CPI*Fcert*Fcf*Fpress; %Total cost of materials
    Ccert = Ceng+Cdev+CFT+Ctool; %Total certification cost
    Cpp = 1035.9*Npp*T^(0.8366)*CPI; %Cost of power plant   
    cost_per_aircraft = Ccert/N
end
