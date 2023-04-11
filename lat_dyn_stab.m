function [Lat_der, Lat_dim_deriv, Lat_modes, A_lat] = lat_dyn_stab(a_El, b_El, bw, sweep, A, Sf, ...
    Sw, Vf, dihedral_angle, CLw, l_f, cw_root, VB, cg, ...
    c_root_tail, bv_tail, Lambda_T, wAC, cw_MAC, theta_tip, M, V_c, AR_T, Sh_tail, c_MAC_tail, CL_tail, bh_tail, x_w, AOA, Inertia, rho, MTOW)
%----------INPUTS-----------
% a_el           = lentgh of long semi-dimension of elliptic fuselage
% b_el           = lentgh of short semi-dimension of elliptic fuselage
% bw             = wing span
% sweep          = sweep angle at c/4 of the wings
% A              = aspect ratio of the wings
% Sf             = surface of the fin
% Sw             = Surface of the wings
% Vf             = Volume fuselage
% dihedral_angle = angle dihedral of the tail
% CLw            = lift coefficient of the wings
% l_f            = length of the fuselage
% cw_root        = root chord of the wing
% VB             = volume of the fuselage
% cg             = position of the center of gravity from the nose
% c_root_tail    = root chord of the tail
% Lambda_T       = sweep angle of the tail, in DEGREES
% wAC            = wing aerodynamic center from the nose
% cw_MAC         = mean aerodynamic chord of the wing
% theta_tip      = twist angle of the wing in RADIANS
% M              = MACH NUMBER
% V_c            = cruise velocity
% AR_T           = aspect ratio of the tail
% Sh_tail        = horizontal surface of the tail
% c_MAC_tail     = mean aerodynamic chord of the tail
% CL_tail        = lift coefficient of the tail
% bh_tail        = horizontal span of the tail
% x_w            = wing position from nose
% AOA            = angle of attack in DEGRESS
show_prints = 0;
%Hypothesis 
% - no wing dihedral angle
% - no delfection angle for flaps
% - cruise conditions with AOA of 2.5???
% - sidewash angle of 2???
alpha = AOA*pi/180; % assuming an angle of attack of 2.5 degrees
alphaF = alpha; % assuming AOA of fuselage of alpha

b_el = 2*b_El;
a_el = 2*a_El;
%% Lateral Stability Derivatives
% Cn_beta -> p.1858, p.1626
bbeta = sqrt(1-M^2);
c1 = 1.1/10 * sin(dihedral_angle); %[deg^-1] see http://airfoiltools.com/airfoil/details?airfoil=naca0010-il
% ds_db = -0.276 + 3.06*Sf/Sw*1/(1+cos(sweep)) + 0.4*Zw/d + 0.0009*A;
SBs = b_el*l_f*0.7; %p.1626, side area of the body (? revoir !!)
if show_prints
    fprintf('Param for p.1633 : xm_lB = %.2f, lB2_SBs = %.2f, sqrt_h1_h2 = %.2f and h_w = %.2f.\n',...
        cg/l_f, l_f^2/SBs, 1, b_el/a_el); %(? revoir !!)
end
KN = 0.0007; %p.1633
nu = 0.000032436;
Re_fus = V_c * l_f/nu;
if show_prints
    fprintf('Reynolds number for fuselage for graph p.1634 is %.5f.\n',Re_fus);
end
KRl = 1.52; %p.1634
Cn_beta_WB = -KN*KRl * SBs/Sw * l_f/bw * 180/pi; %p.1626, per RADIANS

sidewash_and_dyn_press_param = 0.724 + 3.06*(Sf/Sw)/(1 + cos(sweep)) + 0.4*(-1/2) + 0.009*A; %this param = (1+ds_db)*qv/q_inf, see p.1754
t_c = 0.14; %see code wing.m, thickness to chord ratio
clalpha_th = 6.99; %p.478, per RADIANS
clalpha_clalphath = 0.97; %p.478
clalpha = 1.05/bbeta * clalpha_clalphath *clalpha_th; % p.472, per RADIANS
k = 0.75; %p.1668, factor bv/2rs evaluated with bv = fin span and 2rs = fuselage depth in region of the v_tail (? revoir !!)

clalphaM = clalpha/bbeta; %p.471 & p.503, per RADIANS
kappa = clalphaM/(2*pi/bbeta); %p.503
param_p549 = A/kappa*(bbeta + (tan(sweep))^2)^(1/2); %p.549
if show_prints
    fprintf('Param p.549 = %.2f.\n',param_p549); % = 2.07 for last check (? revoir !!)
end
CL_alpha_v = A*1; %p.503 + p.549, per RADIANS
delta_Cy_beta_VWBH = -k*(CL_alpha_v)*sidewash_and_dyn_press_param*Sf/Sw; %p.1645, per RADIANS
zp = bv_tail/4;
lp = l_f - cg - 3/4*c_root_tail + zp/tand(Lambda_T); %see graph p.2777
Cn_beta = Cn_beta_WB - delta_Cy_beta_VWBH*lp/bw; %p.1858, per RADIANS

%% Cl_beta -> p.1845, p.1602
Clbeta_CL_Lambda = -0.001; %p.1563
KM_Lambda = 1.2; %p.1564
lf_param = x_w + sin(sweep)*bw/2; %p.1625
b_param = cos(sweep)*bw; %p.1625
if show_prints
    fprintf('Param p.1625 lf_b = %.2f.\n', lf_param/b_param);
end
Kf = 0.95; %p.1625
Clbeta_CL_A = 0.0001; %Figure p.1564
zW = - b_el/2; %p.1603
d = sqrt(a_el/2 * b_el/2 * pi/0.7854); %p.1603
delta_Clbeta_zW = 1.2*sqrt(A)/57.3 * (zW/b_param)*(2*d/b_param); %p.1603
delta_Clbeta_theta_tansweep = -0.0000325; %see figure p.1566
Cl_beta_WB = (CLw *(Clbeta_CL_Lambda * KM_Lambda * Kf + Clbeta_CL_A) + 0 + delta_Clbeta_zW + theta_tip*180/pi*tan(sweep)*delta_Clbeta_theta_tansweep)*180/pi; %p.1602, per RADIANS
Cl_beta = Cl_beta_WB + delta_Cy_beta_VWBH*(zp*cos(alpha) - lp*sin(alpha))/bw; %p.1845

%% Cy_beta -> p.1812, p.1582
Munk = 0.92; %p.834
BRA = b_el/2 * a_el/2 * pi; %Body Reference Area
CL_alpha_B = 2*Munk*BRA/VB^(2/3); %p.815 of DATCOM
Cy_beta_B = - CL_alpha_B;
delta_Cy_beta_Gam = 0; %p.1582, no dihedral for wings
Ki = -1; %p.1588
Cy_beta_WB = Ki*Cy_beta_B*(BRA/Sw) + delta_Cy_beta_Gam; %p.1582
Cy_beta = Cy_beta_WB + delta_Cy_beta_VWBH; %p.1812

%% Cl_p

%--------- WING-------------
Lambda_beta = atan(tan(sweep)/bbeta);%p.2532
if show_prints
    fprintf('WING : Param Lambda_beta p.2532 = %.3f and the other is %.3f.\n', Lambda_beta*180/pi, bbeta*A/kappa); %last check, 20.566 and 3.151 (? revoir !!)
end
bClp_kappa = -0.235; %p.2550, with param from the "frpintf" just before, per RADIANS
% CLalpha_CL = CL_alpha_v%p.566 same as p.549
CLalpha_CL0 = CL_alpha_v; %here assuming that there's no difference for CL_alpha with CL, per RADIANS
Clp_Cdl_CL2 = -0.01; %p.2554
L = 1.2; %p.723 bcse t_c_max is at 0.37 of chord (see data from http://airfoiltools.com/airfoil/details?airfoil=sc20714-il)
S_wet = Sw - cw_root*a_el;
R_LS = 1.2; %p.749

Re_ell = cw_MAC*V_c/nu;
if show_prints
    fprintf('WING : Reynolds at 30000ft is %.5f.\n',Re_ell); %last check 7e6
end
Cf = 0.003; %p.747
CD0 = Cf*(1 + L*t_c + 100*t_c^4)*R_LS*S_wet/Sw; %p.723
delta_Clp_drag = Clp_Cdl_CL2 * CLw^2 - 1/8 * CD0; %p.2533, per RADIANS
Cl_p_WB = bClp_kappa * kappa/bbeta * CL_alpha_v/CLalpha_CL0 * 1 + delta_Clp_drag; %p.2532 (factor 1 is =1 bcse no dihedral for wings, see eq.7.1.2.2-b), per RADIANS
%---------HORIZONTAL TAIL------------
t_c_tail = 0.10;
Lambda_beta_tail = atan(tand(Lambda_T)/bbeta);%p.2532
if show_prints
    fprintf('TAIL : Param Lambda_beta p.2532 = %.3f and the other is %.3f.\n', Lambda_beta_tail*180/pi, bbeta*AR_T/kappa); %last check, 20.566 and 4.527 (? revoir !!)
end
bClp_kappa_tail = -0.15; %p.2550, with param from the "frpintf" just before
% CLalpha_CL = CL_alpha_v%p.566 same as p.549
CLalpha_CL0_tail = c1; %here assuming that there's no difference for CL_alpha with CL
Clp_Cdl_CL2_tail = -0.02; %p.2554
L = 1.2; %p.723 bcse t_c_max is at 0.30 of chord (see data from http://airfoiltools.com/airfoil/details?airfoil=naca0010-il)
S_wet_tail = Sh_tail;
R_LS = 1.2; %p.749
nu = 0.000032436;
Re_ell_tail = c_MAC_tail*V_c/nu;
if show_prints
    fprintf('TAIL : Reynolds at 30000ft is %.5f.\n',Re_ell_tail); %last check 5.2e6
end
Cf_tail = 0.00325; %p.747
CD0_tail = Cf_tail*(1 + L*t_c_tail + 100*t_c_tail^4)*R_LS*S_wet_tail/Sh_tail; %p.723
delta_Clp_drag_tail = Clp_Cdl_CL2_tail * CL_tail^2 - 1/8 * CD0_tail; %p.2533
z = (zp - lp*tan(alpha))*cos(alpha); %see graph p.2777
clpgamma_clpgamma0_tail = 1 - z/(b_param/2)*sin(dihedral_angle) + 3*(z/(b_param/2))^2*(sin(dihedral_angle))^2;
Cl_p_tail = bClp_kappa_tail * kappa/bbeta * c1/CLalpha_CL0_tail * clpgamma_clpgamma0_tail + delta_Clp_drag_tail; %p.2532 (factor 1 is =1 bcse no dihedral for wings, see eq.7.1.2.2-b), per RADIANS
%--------- Cl_p TOTAL----------------
Cl_p = Cl_p_WB + 0.5*Cl_p_tail*(Sh_tail/Sw)*(bh_tail/bw)^2 + (2*z/bw * (z - zp)/bw)*delta_Cy_beta_VWBH; %p.2783, per RADIANS

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
Cyp_CL_CL0_M0 = 0.05; %see figure p.2529
Cyp_CL_CL0_M = (A + 4*cos(sweep))/(A*B + 4*cos(sweep)) * (A*B + cos(sweep))/(A + cos(sweep)) * Cyp_CL_CL0_M0; %p.2522
Cy_p_WB = K*(Cyp_CL_CL0_M*CLw) + delta_Cyp_gamma; %p.2522
Cy_p = Cy_p_WB + 2*(z - zp)/bw * delta_Cy_beta_VWBH; %p.2776

%% Cn_r -> p.2805, p.2715, p.2593
if show_prints
    fprintf('Param x_bar/c_bar for graph p.2597 is %.2f.\n', x_bar/cw_MAC);
end
Cnr_CL2 = -0.01; %see Figure p.2597
Cnr_CD0 = -0.37; %p.2598
Cn_r_WB = Cnr_CL2 * CLw^2 + Cnr_CD0 * CD0; %p.2593, per RADIANS
Cn_r = Cn_r_WB + 2/bw^2 * (lp*cos(alpha) + zp*sin(alpha))^2 * delta_Cy_beta_VWBH; %p.2805

%% Cl_r -> p.2802, p.2714, p.2581
Clr_CL_CL0_M0 = 0.26; %see Figure p.2589
Clr_CL_CL0_M = (1 + (A*(1-B^2))/(2*B*(A*B + 2*cos(sweep))) + (A*B + 2*cos(sweep))/(A*B + 4*cos(sweep))*((tan(sweep))^2)/8)/(1 + (A + 2*cos(sweep))/(A + 4*cos(sweep))*((tan(sweep))^2)/8)*Clr_CL_CL0_M0; %p.2581
Clbeta_CL = Cl_beta / CLw; %p.1538
delta_Clr_CL = CLw*Clbeta_CL - 0; %p.2582
delta_Clr_theta = 0.0034; %see Figure p.2590
Cl_r_WB = CLw * Clr_CL_CL0_M + delta_Clr_CL + 0 + delta_Clr_theta * theta_tip*180/pi + 0; %p.2581, assuming no dihedral for wings and no delfection of flaps
Cl_r = Cl_r_WB - 2/bw^2 * (lp*cos(alpha) + zp*sin(alpha))*(zp*cos(alpha) - lp*sin(alpha))*delta_Cy_beta_VWBH; %p.2802

%% Cy_r -> p.2799, p.2578
Cy_r_WB = 0; %p.2459 (very bad approx), neglected
Cy_r = Cy_r_WB - 2/bw * (lp*cos(alpha) + zp*sin(alpha))*delta_Cy_beta_VWBH; %p.2799

%% Cy_beta_dot -> p.2825
if show_prints
    fprintf('Param for p.2830 zV_b is %.3f and the one for sigma_beta_WB is %.3f.\n', (zp*cos(alphaF) - lp*(sin(alphaF)))/(b_param/2), (a_el/b_el)/2/(b_param/2));
end
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
% Cl_zeta , see slide 66 course 07
%% Cy_zeta -> slide 64 course 07
% K_prime = %p.1910
% neta = %p.2029
% Kb = %slide 65
% Cy_zeta = c1*alphadeltacL_alphadeltacl * alphadeltacl*K_prime * Kb * Sf/Sw; %slide 64

% Table for lateral stability derivatives
Lat_der = table(Cn_beta, Cl_beta, Cy_beta, Cn_p, Cl_p, Cy_p, Cn_r, Cl_r, Cy_r, Cy_beta_dot, Cl_beta_dot, Cn_beta_dot,'RowNames',{'For AOA = 2.5 deg and cruise cdtn'});

% Computation of A matrix
V0 = V_c;
Ix = Inertia.Ix;
Iy = Inertia.Iy;
Iz = Inertia.Iz;
Ixz = Inertia.Ixz;

Yv = Cy_beta *(1/2*rho*V0*Sw);
Yp = Cy_p *(1/2*rho*V0*Sw*bw);
Yr = Cy_r *(1/2*rho*V0*Sw*bw);
% Yksi=-0.0159*(1/2*rho*V0^2*Sw);
% Yzeta=0.1193*(1/2*rho*V0^2*Sw);
Lv = Cl_beta *(1/2*rho*V0*Sw*bw);
Lp = Cl_p *(1/2*rho*V0*Sw*bw^2);
Lr = Cl_r *(1/2*rho*V0*Sw*bw^2);
% Lksi=0.0454*(1/2*rho*V0^2*Sw*bw);
% Lzeta=0.0086*(1/2*rho*V0^2*Sw*bw);
Nv = Cn_beta *(1/2*rho*V0*Sw*bw);
Np = Cn_p *(1/2*rho*V0*Sw*bw^2);
Nr = Cn_r *(1/2*rho*V0*Sw*bw^2);
% Nksi=0.00084*(1/2*rho*V0^2*Sw*bw);
% Nzeta=-0.0741*(1/2*rho*V0^2*Sw*bw);
Lat_dim_deriv = table(Yv, Yp, Yr, Lv, Lp, Lr, Nv, Np, Nr);
gamae = 0;
thetae = AOA*pi/180 + gamae;
We = V0*sin(thetae);
Ue = V0*cos(thetae);
g = 9.81; %[m/s^2]

M_lat = [MTOW 0 0 0 0;0 Ix -Ixz 0 0;0 -Ixz Iz 0 0;0 0 0 1 0;0 0 0 0 1];
K = [-Yv -(Yp+MTOW*We) -(Yr-MTOW*Ue) -MTOW*g*cos(thetae) -MTOW*g*sin(thetae);
        -Lv -Lp -Lr 0 0; -Nv -Np -Nr 0 0;0 -1 0 0 0;0 0 -1 0 0];
% F=[Yksi Yzeta;Lksi Lzeta;Nksi Nzeta;0 0;0 0];

A_lat = -M_lat\K;
% B_lat = M_lat\F;

eig_lat = eig(A_lat);
disp('eigen values of A_lat :');
disp(eig_lat);
if real(round(eig_lat,10)) <= 0
    LAT_STAB = 'OK';
else 
    LAT_STAB = 'NOT OK';
end

fprintf('Lateral stability is %s\n',LAT_STAB);

%% Lateral Modes of vibrations, according to M.V. COOK
% Spiral mode -> p.216
Ts = -V_c*(Cl_beta*Cn_p - Cl_p*Cn_beta)/(g*(Cl_r*Cn_beta - Cl_beta*Cn_r));
% Roll subsidence -> p.214
Tr = -(Ix*Iz - Ixz^2)/(Iz*Lp + Ixz*Np);
% Dutch roll -> p.217
omega_d = sqrt(Nr*Yv/(Iz*MTOW) + V_c*Nv/Iz);
damp_ratio = - (Nr/Iz + Yv/MTOW)*1/(2*omega_d);
Lat_modes = table(Ts, Tr, omega_d, damp_ratio);