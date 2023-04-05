function [Cn_beta, Cl_beta, Cy_beta, Cn_p, Cl_p, Cy_p, Cn_r, Cl_r, Cy_r, Cy_beta_dot, Cl_beta_dot, Cn_beta_dot] = lat_dyn_stab(a_el, b_el, bw, sweep, A, Sf, ...
    Sw, Vf, dihedral_angle, CLw, l_f, cw_root, VB, cg, ...
    c_root_tail, bv_tail, Lambda_T, wAC, cw_MAC, theta_tip, M, V_c, AR_T, Sh_tail, c_MAC_tail, CL_tail, bh_tail, x_w)
% a_el  = lentgh of long dimension of elliptic fuselage
% b_el  = lentgh of short dimension of elliptic fuselage
% bw    = wing span
% sweep = sweep angle at c/4 of the wings
% A     = aspect ratio of the wings
% Sf    = surface of the fin
% Sw    = Surface of the wings

%Hypothesis 
% - no wing dihedral angle
% - no delfection angle for flaps
% - cruise conditions with AOA of 2.5°
% - sidewash angle of 2°
alpha = 2.5*pi/180; % assuming an angle of attack of 2.5°
beta = 2*pi/180; % assuming a sidewash angle of 2° (beta < 4°)
alphaF = alpha; % assuming AOA of fuselage of alpha

%% Lateral Stability Derivatives
% Cn_beta -> p.1858, p.1626
c1 = 1.1/10 * sin(dihedral_angle); %[deg^-1] see http://airfoiltools.com/airfoil/details?airfoil=naca0010-il
% ds_db = -0.276 + 3.06*Sf/Sw*1/(1+cos(sweep)) + 0.4*Zw/d + 0.0009*A;
SBs = a_el*l_f*0.7; %p.1626, side area of the body (à revoir !!)
fprintf("Param for p.1633 : xm_lB = %.2f, lB2_SBs = %.2f, sqrt_h1_h2 = %.2f and h_w = %.2f.\n",...
    cg/l_f, l_f^2/SBs, 1, a_el/b_el); %(à revoir !!)
KN = 0.0009; %p.1633
nu = 0.000032436;
Re_fus = V_c * l_f/nu;
fprintf("Reynolds number for fuselage for graph p.1634 is %.5f.\n",Re_fus);
KRl = 1.52; %p.1634
Cn_beta_WB = -KN*KRl * SBs/Sw * l_f/bw * 180/pi; %p.1626, per RADIANS

sidewash_and_dyn_press_param = 0.724 + 3.06*(Sf/Sw)/(1 + cos(sweep)) + 0.4*(-1/2) + 0.009*A; %this param = (1+ds_db)*qv/q_inf, see p.1754
t_c = 0.14; %see code wing.m, thickness to chord ratio
clalpha = 6.28 + 4.7*t_c; % p.471
clalphaM = clalpha/beta;
kappa = clalphaM/(2*pi/beta);
param_p549 = A/kappa*(beta + (tan(sweep))^2)^(1/2);
fprintf("Param p.549 = %.2f.\n",param_p549); % = 2.07 for last check (à revoir !!)
k = 0.75; %p.1668, factor bv/2rs evaluated with bv = fin span and 2rs = fuselage depth in region of the v_tail (à revoir !!)
CL_alpha_v = A*1.3;%p.503 + p.549
delta_Cy_beta_VWBH = -k*(CL_alpha_v)*sidewash_and_dyn_press_param*Sf/Sw; %p.1645, per RADIANS
zp = bv_tail/4;
lp = l_f - cg - 3/4*c_root_tail + zp/tand(Lambda_T); %see graph p.2777
Cn_beta = Cn_beta_WB + delta_Cy_beta_VWBH*lp/bw; %p.1858, per RADIANS

%% Cl_beta -> p.1845, p.1602
Clbeta_CL_Lambda = -0.001; %p.1563
KM_Lambda = 1.22; %p.1564
lf_param = x_w + sin(sweep)*bw/2;
b_param = cos(sweep)*bw;
fprintf("Param p.1625 lf_b = %.2f.\n", lf_param/b_param);
Kf = 0.94; %p.1625
Clbeta_CL_A = 0.0001; %p.1564
zW = - a_el/2; %p.1603
d = sqrt(a_el/2 * b_el/2 * pi/0.7854);
delta_Clbeta_zW = 1.2*sqrt(A)/57.3 * (zW/b_param)*(2*d/b_param); %p.1603
delta_Clbeta_theta_tansweep = 0.0000325; %see figure p.1566
Cl_beta_WB = (CLw *(Clbeta_CL_Lambda * KM_Lambda * Kf + Clbeta_CL_A) + 0 + delta_Clbeta_zW + theta_tip*180/pi*tan(sweep)*delta_Clbeta_theta_tansweep)*180/pi; %p.1602, per RADIANS
Cl_beta = Cl_beta_WB + delta_Cy_beta_VWBH*(zp*cos(alpha) - lp*sin(alpha))/bw; %p.1845

%% Cy_beta -> p.1512, p.1582
Munk = 0.95; %p.834
BRA = b_el/2 * a_el/2 * pi; %Body Reference Area
CL_alpha_B = 2*Munk*BRA/VB^(2/3); %p.815 of DATCOM
Cy_beta_B = - CL_alpha_B;
delta_Cy_beta_Gam = 0; %p.1582, no dihedral for wings
Ki = -1; %p.1588
Cy_beta_WB = Ki*Cy_beta_B*(BRA/Sw) + delta_Cy_beta_Gam; %p.1582
Cy_beta = Cy_beta_WB + delta_Cy_beta_VWBH; %p.1812

%% Cl_p
bbeta = sqrt(1-M^2);
%--------- WING-------------
Lambda_beta = atan(tan(sweep)/bbeta);%p.2532
fprintf("WING : Param Lambda_beta p.2532 = %.3f and the other is %.3f.\n", Lambda_beta*180/pi, bbeta*A/kappa); %last check, 20.566 and 4.527 (à revoir !!)
bClp_kappa = -0.28; %p.2550, with param from the "frpintf" just before
% CLalpha_CL = CL_alpha_v%p.566 same as p.549
CLalpha_CL0 = CL_alpha_v; %here assuming that there's no difference for CL_alpha with CL
Clp_Cdl_CL2 = -0.012; %p.2554
L = 1.2; %p.723 bcse t_c_max is at 0.37 of chord (see data from http://airfoiltools.com/airfoil/details?airfoil=sc20714-il)
S_wet = Sw - cw_root*b_el;
R_LS = 1.2; %p.749

Re_ell = cw_MAC*V_c/nu;
fprintf("WING : Reynolds at 30000ft is %.5f.\n",Re_ell); %last check 7e6
Cf = 0.0031; %p.747
CD0 = Cf*(1 + L*t_c + 100*t_c^4)*R_LS*S_wet/Sw; %p.723
delta_Clp_drag = Clp_Cdl_CL2 * CLw^2 - 1/8 * CD0; %p.2533
Cl_p_WB = bClp_kappa * kappa/bbeta * CL_alpha_v/CLalpha_CL0 * 1 + delta_Clp_drag; %p.2532 (factor 1 is =1 bcse no dihedral for wings, see eq.7.1.2.2-b)
%---------HORIZONTAL TAIL------------
t_c_tail = 0.10;
Lambda_beta_tail = atan(tand(Lambda_T)/bbeta);%p.2532
fprintf("TAIL : Param Lambda_beta p.2532 = %.3f and the other is %.3f.\n", Lambda_beta_tail*180/pi, bbeta*AR_T/kappa); %last check, 20.566 and 4.527 (à revoir !!)
bClp_kappa_tail = -0.2; %p.2550, with param from the "frpintf" just before
% CLalpha_CL = CL_alpha_v%p.566 same as p.549
CLalpha_CL0_tail = c1; %here assuming that there's no difference for CL_alpha with CL
Clp_Cdl_CL2_tail = -0.02; %p.2554
L = 1.2; %p.723 bcse t_c_max is at 0.30 of chord (see data from http://airfoiltools.com/airfoil/details?airfoil=naca0010-il)
S_wet_tail = Sh_tail;
R_LS = 1.2; %p.749
nu = 0.000032436;
Re_ell_tail = c_MAC_tail*V_c/nu;
fprintf("TAIL : Reynolds at 30000ft is %.5f.\n",Re_ell_tail); %last check 5e6
Cf_tail = 0.00325; %p.747
CD0_tail = Cf_tail*(1 + L*t_c_tail + 100*t_c_tail^4)*R_LS*S_wet_tail/Sh_tail; %p.723
delta_Clp_drag_tail = Clp_Cdl_CL2_tail * CL_tail^2 - 1/8 * CD0_tail; %p.2533
Cl_p_tail = bClp_kappa_tail * kappa/bbeta * c1/CLalpha_CL0_tail * 1 + delta_Clp_drag_tail; %p.2532 (factor 1 is =1 bcse no dihedral for wings, see eq.7.1.2.2-b)
%--------- Cl_p TOTAL----------------
z = (zp - lp*tan(alpha))*cos(alpha); %see graph p.2777
Cl_p = Cl_p_WB + 0.5*Cl_p_tail*(Sh_tail/Sw)*(bh_tail/bw)^2 + (2*z/bw * (z - zp)/bw)*delta_Cy_beta_VWBH; %p.2783

%% Cn_p -> p.2793, p.2711, p.2559, p.2531 (for Cl_p)
x_bar = wAC - cg;%p.2560
B = sqrt(1 - M^2*(cos(sweep))^2);%p.2559
K = 0.1; %very arbitrary !! see p.2524 but too much caclulations !
Cnp_CL_CL0_M0 = - 1/6 * (A + 6*(A + cos(sweep))*(x_bar/cw_MAC * tan(sweep)/A + ((tan(sweep))^2)/12))/(A + 4*cos(sweep)); %p.2559
Cnp_CL_CL0_M = (A + 4*cos(sweep))/(A*B + 4*cos(sweep))*((A*B + 1/2*(A*B + cos(sweep))*(tan(sweep))^2)/(A + 1/2*(A + cos(sweep))*(tan(sweep))^2))*Cnp_CL_CL0_M0; %p.2559
delta_Cnp_theta = 0.00025; %see Figure p.2569
Cn_p_WB = - Cl_p*tan(alpha) - K*(-Cl_p*tan(alpha) - Cnp_CL_CL0_M * CLw) + delta_Cnp_theta*theta_tip*180/pi; %simplified because cruise condition => no deflection angle for flaps (delta_f), p.2559
Cn_p = Cn_p_WB - 2/bw * (lp*cos(alpha) + zp*sin(alpha))*(z-zp)/bw * delta_Cy_beta_VWBH;%p.2793

%% Cy_p -> p.2776, p.2696, p.2522
delta_Cyp_gamma = 0; %p.2522 bcse no dihedral angle for wings
Cyp_CL_CL0_M0 = 0.07; %see figure p.2529
Cyp_CL_CL0_M = (A + 4*cos(sweep))/(A*B + 4*cos(sweep)) * (A*B + cos(sweep))/(A + cos(sweep)) * Cyp_CL_CL0_M0; %p.2522
Cy_p_WB = K*(Cyp_CL_CL0_M*CLw) + delta_Cyp_gamma; %p.2522
Cy_p = Cy_p_WB + 2*z/bw * (z - zp)/bw * delta_Cy_beta_VWBH; %p.2776

%% Cn_r -> p.2805, p.2715, p.2593
fprintf("Param x_bar/c_bar for graph p.2597 is %.2f.\n", x_bar/cw_MAC);
Cnr_CL2 = 0.01; %see Figure p.2597
Cnr_CD0 = -0.33; %p.2598
Cn_r_WB = Cnr_CL2 * CLw^2 + Cnr_CD0 * CD0; %p.2593
Cn_r = Cn_r_WB + 2/bw^2 * (lp*cos(alpha) + zp*sin(alpha))^2 * delta_Cy_beta_VWBH; %p.2805

%% Cl_r -> p.2802, p.2714, p.2581
Clr_CL_CL0_M0 = 0.25; %see Figure p.2589
Clr_CL_CL0_M = (1 + (A*(1-B^2))/(2*B*(A*B + 2*cos(sweep))) + (A*B + 2*cos(sweep))/(A*B + 4*cos(sweep))*((tan(sweep))^2)/8)/(1 + (A + 2*cos(sweep))/(A + 4*cos(sweep))*((tan(sweep))^2)/8)*Clr_CL_CL0_M0; %p.2581
Clbeta_CL = Cl_beta / CLw;%p.1538
delta_Clr_CL = CLw*Clbeta_CL - 0; %p.2582
delta_Clr_theta = 0.0034; %see Figure p.2590
Cl_r_WB = CLw * Clr_CL_CL0_M + delta_Clr_CL + 0 + delta_Clr_theta * theta_tip*180/pi + 0; %p.2581, assuming no dihedral for wings and no delfection of flaps
Cl_r = Cl_r_WB - 2/bw^2 * (lp*cos(alpha) + zp*sin(alpha))*(zp*cos(alpha) - lp*sin(alpha))*delta_Cy_beta_VWBH; %p.2802

%% Cy_r -> p.2799, p.2578
Cy_r_WB = 2.6*0; %p.2459 (very bad approx), neglected
Cy_r = Cy_r_WB - 2/bw * (lp*cos(alpha) + zp*sin(alpha))*delta_Cy_beta_VWBH; %p.2799

%% Cy_beta_dot -> p.2825
fprintf("Param for p.2830 zV_b is %.3f and the one for sigma_beta_WB is %.3f.\n", (zp*cos(alphaF) - lp*(sin(alphaF)))/(b_param/2), (a_el/b_el)/2/(b_param/2));
sigma_beta_alpha = (-0.002-0.055)/2; %p.2840 + p.2841, mean value
sigma_beta_theta = (-0.0055 - 0.05)/2; %p.2860 + p.2861, same param of p.2830 (mean value)
sigma_beta_WB = (0.204 - 0.186)/2; %p.2880 + p.2881
sigma_beta = sigma_beta_alpha *alphaF*180/pi + 0 - sigma_beta_theta*theta_tip*180/pi + sigma_beta_WB; %p.2826
Cy_beta_dot = 2*c1 * sigma_beta *Sf/Sw * (lp*cos(alphaF) + zp*sin(alphaF))/bw; %p.2825


%% Cl_beta_dot -> p.2886
Cl_beta_dot = Cy_beta_dot * (zp*cos(alphaF) - lp*sin(alphaF))/bw;

%% Cn_beta_dot -> p.2888
Cn_beta_dot = - Cy_beta_dot * (lp*cos(alphaF) + zp*sin(alphaF))/bw;

% Cn_zeta =
% Cl_zeta =
% Cy_zeta =