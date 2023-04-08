function [T] = long_dyn_stab(MTOW,...
    a_el,b_el,bw,Sw,CLw_alpha,rho,V_c,ARw,M,Altitude,CL_alphaT,Sh_tail,de_dAOA1,...
    static_stability,AOA,alpha_L0,l_f,l_cg,sweep,cl_alphaw,...
    l_arm,V_hT,cw_MAC,Xw,cw_root,xw_AC)
% Parameters
show_prints = 0;
    % CL_alpha
d = (2*a_el+2*b_el)/2;
K_WB = 1-0.25*(d/bw)^2+0.025*(d/bw);
CL_alphaWB = K_WB*CLw_alpha;
heta_H = 0.95; %[0.9;1]

    % CD_alpha
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

% if show_prints
%     fprintf('CD_DM = %.2f.\n',CD_DM);
%     fprintf('CD_M = %.2f.\n',CD_M);
%     fprintf('DM = %.2f.\n',DM);
% end
CD_Cm = (CD_DM+CD_M)/DM;
    % Cm_u
% param1 = tan(sweep)/sqrt(1-M^2);
% param2 = tan(sweepDM)/sqrt(1-MDM^2);
% if show_prints
%     fprintf('Param1 slide 31 = %.2f.\n',param1);
% end
% if show_prints
%     fprintf('Param2 slide 31 = %.2f.\n',param2);
% end

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
M0=0;
beta0 = sqrt(1-M0^2);
k1 = cl_alphaw/(2*pi);
CL0_alpha_wing = 2*pi*ARw/(2+sqrt(((ARw^2*beta0^2)/k1^2)*(1+tan(sweep)^2/beta0^2)+4));
B = sqrt(1-M^2*cos(sweep)^2);
CL_qw = (ARw+2*cos(sweep))/(ARw*B+2*cos(sweep))*(1/2+2*Xw/cw_MAC)*CL0_alpha_wing;
CL_qH = 2*CL_alphaT*heta_H*V_hT;

    % Cm_q
K = 0.72; % From Figure 5 in FD&C L.6
Cm_qw0 = -K*cl_alphaw*cos(sweep) * (ARw*(2*(Xw/cw_MAC)^2+1/2*(Xw/cw_MAC))/(ARw+2*cos(sweep)) + 1/24*ARw^3*tan(sweep)^2/(ARw+6*cos(sweep)) + 1/8);
Cm_qw = Cm_qw0 * ( ( (ARw^3*tan(sweep)^2)/(ARw*B+6*cos(sweep)) +3/B) / ((ARw^3*tan(sweep)^2)/(ARw+6*cos(sweep))+3) );
Cm_qH = -2*CL_alphaT*heta_H*V_hT*l_arm/cw_MAC;

    % CL_adot
CLg = -3.88;
CL_adotw = 1.5*xw_AC/cw_root*CLw_alpha + 3*CLg;
CL_adotH = 2*CL_alphaT*heta_H*V_hT*de_dAOA1;

    % Cm_adot
Cm_adotw = 0;
Cm_adotH = -2*CL_alphaT*heta_H*V_hT*l_arm/cw_MAC*de_dAOA1;

    % CL_df
beta = sqrt(1-M^2);
alphaCL_cl = 1.05; %cf/c = 0.4
cld_th = 5; %approximated don't know whats t/c [1/rad]
cl_cld_th = 0.9; %approximated wrong value of cl_alpha
K3 = 1; % = 1 if df<12 [deg]
cl_df = cld_th*cl_cld_th*K3;
Kb = 0.33; %very approximated

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

% fprintf('CL_qw = %.4d and CL_qH = %.4d\n',CL_qw, CL_qH);
CL_q = CL_qw + CL_qH;
CD_q = 0; %usually neglected
Cm_q = Cm_qw + Cm_qH;

CL_adot = CL_adotw + CL_adotH;
CD_adot = 0; %usually neglegted
Cm_adot = Cm_adotw + Cm_adotH;

CL_df = cl_df*beta*CLw_alpha/cl_alphaw*alphaCL_cl*Kb;
Cm_df = 0; %neglected

CL_ih = CL_alphaT*Sh_tail/Sw;
Cm_ih = -CL_alphaT*V_hT;

CL_heta = CL_df*Sh_tail/Sw;
Cm_heta = -CL_df*V_hT;

T = table(CL_alpha, CD_alpha, Cm_alpha, CL_u, CD_u, Cm_u, CL_q, CD_q, Cm_q, CL_adot, CD_adot, Cm_adot, CL_df, Cm_df, CL_ih, Cm_ih, CL_heta, Cm_heta,'RowNames',{'For AOA = 2.5 deg and cruise cdtn'});


end

