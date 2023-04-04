function [CL_alpha,CD_alpha,Cm_alpha,CL_u,CD_u,Cm_u,CL_q,CD_q,Cm_q,CL_adot,CD_adot,Cm_adot,CL_df,Cm_df,CL_ih,Cm_ih,CL_heta,Cm_heta] = long_dyn_stab(MTOW,a_el,b_el,bw,Sw,CLw_alpha,rho,V_c,ARw,M,Altitude,CL_alphaT,Sh_tail,de_dAOA1,static_stability,AOA,alpha_L0,l_f,l_cg)
% Parameters

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
[bwDM,SwDM,CLw_alphaDM,CDw_alphaDM,CLwDM,CDwDM,DDM,cw_rootDM,cw_tipDM,cw_MACDM,xw_ACDM,yw_ACDM,Vw_fuelDM,sweepDM,cDM,alpha_L0DM] = wing(MDM,Altitude,0.95*MTOW,AOA);
[S_tailDM,Sh_tailDM,Sv_tailDM,c_root_tailDM,c_tip_tailDM, dihedral_angleDM, l_armDM, CL_tailDM, Lambda_TDM,...
    b_tailDM, bv_tailDM, bh_tailDM, W_tailDM, Cn_beta_AhDM, V_vfDM, hight_rootDM, hight_tipDM,...
    rudder_chord_rootDM, rudder_chord_tipDM, rudder_chordDM, S_rudderDM] = v_tail(MTOW,...
    2*b_el,cw_MACDM,sweepDM*180/pi,SwDM,l_f,l_cg,bwDM);
[W_wingDM, W_fuselageDM, W_landing_gear_noseDM, W_landing_gear_mainDM, W_installed_engineDM, W_payloadDM, W_FSDM, W_fuelDM, W_systemDM, W_totDM, W_subsystDM, W_sensorsDM] = mass(MTOW,bwDM,cw_rootDM,cw_tipDM,l_armDM);
WDM = [W_fuselageDM;W_wingDM;W_tailDM;W_installed_engineDM;W_landing_gear_noseDM;W_landing_gear_mainDM;W_payloadDM;W_fuelDM+W_FSDM;W_subsystDM;W_sensorsDM];
MTOWDM = sum(WDM);
CL = MTOWDM*9.81/(1/2*rhoDM*V_cDM^2*SwDM);
CD_DM = 0.017+CL^2/(0.8*pi*ARw);
CD_Cm = (CD_DM+CD_M)/DM;

    % CM_u
x_ACwing = 0;
dxac_dm = (xw_ACDM-x_ACwing)/DM;

    % CL_q
A = 0;
c = 0;
CLaw = 0; %with M=0
Lambda = 0;
Xw = 0;
B = 0;
CL_qw = (A+2*cos(Lambda))/(A*B+2*cos(Lambda))*(1/2+2*Xw/c)*CLw_alpha;
heta = 0;
V_T = 0;
CL_qH = 2*CL_alphaT*heta*V_T;

%Longitudinal Derivatives (for the entire aircraft) [rad^-1]
CL_alpha = CL_alphaWB + CL_alphaT*heta_H*Sh_tail/Sw*(1-de_dAOA1);
CD_alpha = 2*CL*CL_alpha/(pi*ARw*e);
Cm_alpha = static_stability;
% CM_alpha = dCM_dCL*CL_alpha;
CL = CL_alpha*(AOA*pi/180-alpha_L0);
CL_u = M^2/(1-M^2)*CL;
CD_u = M*CD_Cm;
Cm_u = -CL+dxac_dm;

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

