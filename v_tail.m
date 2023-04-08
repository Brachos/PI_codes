function [S_tail, S_h,S_v,c_root_tail, c_tip_tail, angle,l_arm,C_L,Lambda_T, ...
    b_tail, b_v, b_h, W_tail,CN_tot, V_v, hight_root, hight_tip, ...
    rudder_chord_root, rudder_chord_tip, rudder_chord, S_rudder, AR_t] = ...
    v_tail(Mass, h_f_max, c_chord, Lambda_LE, Sw, l_f, l_cg, bw)
%Code destin? ? obtenir les principaux param?tres g?om?trique de la tail en
%fonction des caract?ristiques des ailes. Cette m?thode est bas?e sur
%l'ouvrage de r?f?rence "Aircraft design, A Systems Engineering Approach"
%de Mohammed H Sadraey et plus pr?cis?ment du chapitre 6.

%Cette m?thode est it?rative et s'appuie sur les ?quilibre statiques
%et dynamique.

show_prints = 1;
%% Figures settings
feature('DefaultCharacterSet','UTF8');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultTextFontsize',13);
set(groot, 'defaultAxesFontsize',13);
set(groot, 'defaultLegendFontsize',13);
set(groot, 'defaultLegendLocation','best');
set(0, 'DefaultLineLineWidth', 1.8);
pt = 12;
%% Data required
c_bar = c_chord;
% - lambda
lambda = 0.3; % fixed
% - sweep angle
Lambda_T = Lambda_LE + 5;
% - aspect ratio
AR_t = 4;
% conversions
pound = 2.20462262; % kg to lbs
feet = 3.28; % m to ft
%% Longitudinal trim
% Horizontal tail

% Le coefficient (l/c_bar * S_h/S)=:V_h_bar  est tr?s important et est
% appel? "horizontal tail volume coefficient". Il joue un role tr?s
% important dans la stabilit? longitutinale de l'appareil. Apr?s avoir
% consult? plusieur sources statistiques, ce cefficient sera fix? ? une
% premi?re valeur de 0.6;
V_h_bar = 0.6;
l_arm = l_f - l_cg - 0.75;
% Ainsi, on peut maintenant d?temriner S_h :
S_h = V_h_bar*c_bar*Sw/l_arm;

% Vertical fin, see slide 57 course Noels Preminilary design
K_beta = 0.3*l_cg/l_f + 0.75 *h_f_max/l_f - 0.105;
S_fs = 0.8*l_f*h_f_max;
CN_beta_f = -K_beta*S_fs*l_f/Sw/bw; %à revoir !!
CN_beta_i = -0.017; %because high wings
CN_tot = CN_beta_f + CN_beta_i;
if show_prints
    fprintf('Param slide 57 for vertical tail is %.3f.\n', CN_tot);
end
% voir graphique slide 57
V_v = 0.1; % avec V_v = S_F*l_F/(S*b)
l_F = l_arm; % first guess, distance between cg and fin ac
S_v = V_v*Sw*bw/l_F;

angle = atan(sqrt(S_v/S_h)); % angle de la v_tail
C_L = 0.15*cos(angle);

S_tail = S_h + S_v;
b_tail = sqrt(AR_t*S_tail);% span along the tail (one side)
b_h = sin(angle)*b_tail;
b_v = cos(angle)*b_tail/2;
lambda_t = 0.6;
c_root_tail = 2*S_tail/b_tail/(1+lambda_t);
c_tip_tail = lambda_t * c_root_tail;
%% Weight
[W_mass] = fuel_weight();
W_fuel = W_mass*pound;
W_dg = Mass*pound - 0.45*W_fuel;
N_z = 1.5*3; %Ultimate load factor = 1.5*limit load factor
Sw = Sw*feet^2; %[ft^2]
q = 8.6292e+03*pound/feet^2;
W_tail_h = 0.016*(N_z*W_dg)^(0.414)*q^(0.168)*(S_tail*feet^2*sin(angle))^(0.896)*...
    (100*0.1/cosd(Lambda_LE))^(-0.12)*(AR_t/(cosd(Lambda_T))^2)^(0.043)*lambda_t^(-0.02);

%Ht/Hv ??? = 0.5 because no information about V-tail
Ht = 1;
Hv = 2;
W_tail_v = 0.073*(1+0.2*Ht/Hv)*(N_z*W_dg)^(0.376)*q^(0.122)*(S_tail*feet^2*cos(angle))^0.873...
    *(100*0.1/cos(Lambda_T))^-0.49*(AR_t/(cos(Lambda_T))^2)^0.357*lambda^0.039;
% disp(W_tail_h);
% disp(W_tail_v);
Sw = Sw/feet^2;

W_tail = (W_tail_v + W_tail_h)/pound;
%% control surfaces
%Based on Raymer tab
span_covered = 0.7; %proportion of total span
rudder_chord = 0.35; %proportion of the chord
hight_root = (1-span_covered)/2*b_v;
hight_tip = ((1-span_covered)/2+span_covered)*b_v;
rudder_chord_root = interp1([0 b_v],[c_root_tail c_tip_tail],hight_root)*rudder_chord;
rudder_chord_tip = interp1([0 b_v],[c_root_tail c_tip_tail],hight_tip)*rudder_chord;
S_rudder = (rudder_chord_root + rudder_chord_tip)*(hight_tip-hight_root)/2;
%% Prints
if show_prints
    fprintf('Tail surface : %.2dft?\n',S_tail);
    fprintf('Tail span : %.2dm\n',b_tail);
    fprintf('Tail horizontal span : %.2dm\n',b_h);
    fprintf('Tail vertical span : %.2dm\n',b_v);
    fprintf('Dihedral angle (degrees) : %.2d\n',angle*180/pi);
    fprintf('Surface ratio : %.2d\n',S_h/Sw);
    fprintf('Rudder surface : %.2d\n', S_rudder);
    fprintf('Weight of the horizontal tail : %.2dkg\n', W_tail_h/pound);
    fprintf('Weight of the vertical tail : %.2dkg\n', W_tail_v/pound);
    fprintf('Total weight of the tail : %.2dkg\n', W_tail);
end
