function [Cn_beta, Cl_beta] = lat_dyn_stab(a_el, b_el, bw, sweep, A, Sf, Sw, Cn_beta_Ah, Vf, dihedral_angle, CLw, l_f, cw_root, VB)
% a_el  = lentgh of long dimension of elliptic fuselage
% b_el  = lentgh of short dimension of elliptic fuselage
% bw    = wing span
% sweep = sweep angle of c/4 of the wings
% A     = aspect ratio of the wings

alpha = 3*pi/180; % assuming an angle of attack of 3°
beta = 2*pi/180; % assuming a sidewash angle of 2° (beta < 4°)

%Lateral Derivatives
% Cn_beta
Zw = b_el/2;
c1 = 1.1/10 * sin(dihedral_angle); %[deg^-1] see http://airfoiltools.com/airfoil/details?airfoil=naca0010-il
d = sqrt(a_el/2 * b_el/2 * pi/0.7854);
ds_db = -0.276 + 3.06*Sf/Sw*1/(1+cos(sweep)) + 0.4*Zw/d + 0.0009*A;
Cn_beta_tail = Vf * c1 * (1 - ds_db);
Cn_beta = Cn_beta_Ah + Cn_beta_tail;
% Cl_beta
Cl_beta_CLw = -0.0005; % slide 51, course lateral stability
Cl_beta = Cl_beta_CLw * CLw;
% Cy_beta 
Munk = 0.95; %p.834
BRA = b_el/2 * a_el/2 * pi; %Body Reference Area
CL_alpha_B = 2*Munk*BRA/VB^(2/3); %p.815 of DATCOM
Cy_beta_B = - CL_alpha_B;
delta_Cy_beta_Gam = 0; %p.1582, no dihedral for wings
Ki = -1; %p.1588 DATCOM
Cy_beta_WB = Ki*Cy_beta_B*(BRA/Sw) + delta_Cy_beta_Gam; %p.1582
% Cn_p = 
% Cl_p =
% Cy_p = 
% Cn_r = 
% Cl_r = 
% Cy_r = 
% Cn_xi =
% Cl_xi =
% Cy_xi =
% Cn_zeta =
% Cl_zeta =
% Cy_zeta =