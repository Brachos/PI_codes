function [CL_alpha,CD_alpha,Cm_alpha,CL_u,CD_u,Cm_u,CL_q,CD_q,Cm_q,CL_adot,...
    CD_adot,Cm_adot,CL_df,Cm_df,CL_ih,Cm_ih,CL_heta,Cm_heta] = long_dyn_stab(MTOW,...
    a_el,b_el,bw,Sw,CLw_alpha,rho,V_c,ARw,M,Altitude,CL_alphaT,Sh_tail,de_dAOA1,...
    static_stability,AOA,alpha_L0,l_f,l_cg,sweep,feet,pound,cl_alphaw,l_arm)
% Parameters
% beta = sqrt(1-M^2);
% k = cl_alpha/(2*pi);
% CL_alpha_wing = 2*pi*ARw/(2+sqrt(((ARw^2*beta^2)/k^2)*(1+tan(sweep)^2/beta^2)+4));
% fprintf('CL_alpha of the wing DATCOM method : %.4d\n',CL_alpha_wing);
% CL_alpa
% CL1 = CL_
    % CL_alpha
d = (2*a_el+2*b_el)/2;
K_WB = 1-0.25*(d/bw)^2+0.025*(d/bw);
CL_alphaWB = K_WB*CLw_alpha;
heta_H = 0.95; %[0.9;1]

    % CD_alpha
V0 = 0;
CL = MTOW*9.81/(1/2*rho*V_c^2*Sw);
CD = 0.017+CL^2/(0.8*pi*ARw);
e = 0.8;

    % CM_alpha
% sum = 0;
% for i=1:Nelem
%     sum = sum+W(i)*de_dAOA1*W^2;
% end
% DM_dalpha = 1/2*rho*V_c^2/36.5*1;
% x_acWB = 0;
% x_ac = (x_acWB+CL_alphaT/CL_alphaWB*heta_H*Sh_tail/S*x_ach*(1-de_dAOA1))/(1+CL_alphaT/CL_alphaWB*heta_H*Sh_tail/Sw*(1-de_dAOA1));
% dCM_dCL = x_cg-x_ac;

    % CD_u
CD_M = CD;
DM = M/20;
MDM = M+DM;
[speedDM,rhoDM] = speed(Altitude,MDM);
V_cDM = speedDM;
[bwDM,SwDM,~,~,~,~,~,cw_rootDM,cw_tipDM,cw_MACDM,~,~,~,sweepDM,~,~] = wing(MDM,...
    Altitude,0.95*MTOW,AOA);

[~,~,~,~,~,~,l_armDM,~,~,~,~,~,W_tailDM,~, ~, ~, ~,~,~, ~, ~] = v_tail(MTOW,...
    2*b_el,cw_MACDM,sweepDM*180/pi,SwDM,l_f,l_cg,bwDM);

[W_wingDM, W_fuselageDM, W_landing_gear_noseDM, W_landing_gear_mainDM,...
    W_installed_engineDM, W_payloadDM, W_FSDM, W_fuelDM, ~, ~,...
    W_subsystDM, W_sensorsDM] = mass(MTOW,bwDM,cw_rootDM,cw_tipDM,l_armDM);

WDM = [W_fuselageDM;W_wingDM;W_tailDM;W_installed_engineDM;W_landing_gear_noseDM;W_landing_gear_mainDM;W_payloadDM;W_fuelDM+W_FSDM;W_subsystDM;W_sensorsDM];
MTOWDM = sum(WDM);
CL_DM = MTOWDM*9.81/(1/2*rhoDM*V_cDM^2*SwDM);
CD_DM = 0.017+CL_DM^2/(0.8*pi*ARw);
fprintf('CD_DM = %.2f.\n',CD_DM);
fprintf('CD_M = %.2f.\n',CD_M);
fprintf('DM = %.2f.\n',DM);
CD_Cm = (CD_DM+CD_M)/DM;

    % CM_u
param1 = tan(sweep)/sqrt(1-M^2);
fprintf('Param1 slide 31 = %.2f.\n',param1);
param2 = tan(sweepDM)/sqrt(1-MDM^2);
fprintf('Param2 slide 31 = %.2f.\n',param2);
K1 = 1.4;
K2 = 0.25;
xac_cr = 0.45;
xac_crDM = 0.45; %They are similar since the slope of the corresponding curve is flat
% Parameters above determined graphically with experimental relations in
% Flight Dynamics and control lesson 6 slide 30
x_ACwing = K1*(xac_cr-K2);
x_ACwingDM = K1*(xac_crDM-K2);
dxac_dm = (x_ACwingDM-x_ACwing)/DM;

    % CL_q
% Computing CL_alpha with M=0
M0=0;
beta = sqrt(1-M0^2);
k = cl_alphaw/(2*pi);
CL0_alpha_wing = 2*pi*ARw/(2+sqrt(((ARw^2*beta^2)/k^2)*(1+tan(sweep)^2/beta^2)+4));
k = cl_alphaT/(2*pi);
[bw0,Sw0,CLw_alpha0,~,~,~,~,~,~,~,~,~,~,~,~,~] = wing(0,Altitude,0.95*MTOW,AOA);
[~,Sh_tail0,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = v_tail(MTOW,2*b_el,...
cw_MACDM,sweepDM*180/pi,SwDM,l_f,l_cg,bwDM);

K_WB = 1-0.25*(d/bw0)^2+0.025*(d/bw0);
CL_alphaWB0 = K_WB*CLw_alpha0;
zwt = 1;
m0 = zwt/(bw0/2);
de_dAOA0 = 1.75*CLw_alpha0/(pi*7*(2*l_arm*0.3/bw0)^(1/4)*(1+m0));
CL0_alpha = CL_alphaWB0 + CL_alphaT*heta_H*Sh_tail0/Sw0*(1-de_dAOA0); %with M=0
fprintf('CL_alpha with M=0 = %.4d\n',CL0_alpha);
c = 0;
Lambda = 0;
Xw = 0;
B = 0;
CL_qw = (ARw+2*cos(Lambda))/(ARw*B+2*cos(Lambda))*(1/2+2*Xw/c)*CL0_alpha;
heta = 0;
V_T = 0;
CL_qH = 2*CL_alphaT*heta*V_T;

%Longitudinal Derivatives (for the entire aircraft) [rad^-1]
CL_alpha = CL_alphaWB + CL_alphaT*heta_H*Sh_tail/Sw*(1-de_dAOA1);
CD_alpha = 2*CL*CL_alpha/(pi*ARw*e);
Cm_alpha = static_stability;
% CM_alpha = dCM_dCL*CL_alpha;
CL = CL_alpha*(AOA*pi/180-alpha_L0); % Correct value
% fprintf('Lift coefficient computed with dynamic stability = %.2f.\n',CL);
% rho_mat = 0.48; %[kg/m^3]
% L = CL*Sw*1/2*V_c^2*rho_mat;
% fprintf('Lift computed with dynamic stability = %.2f.\n',L);
CL_u = M^2/(1-M^2)*CL;
CD_u = M*CD_Cm;
Cm_u = -CL*dxac_dm;

CL_q = CL_qw + CL_qH;
CD_q = 0; %usually neglected
Cm_q = 0;

CL_adot = 0;
CD_adot = 0;
Cm_adot = 0;

CL_df = 0;
Cm_df = 0;

CL_ih = 0;
Cm_ih = 0;

CL_heta = 0;
Cm_heta = 0;


end

