function [S_tail, S_h,S_v,c_root_tail, c_tip_tail, angle,l_arm,C_L,Lambda_T, b_tail, b_v, b_h, W_tail] = v_tail(Mass,...
    D_f_max,h_f_max,V_c,c_chord,Lambda_LE,S,l_f,l_cg,b)
%Code destin? ? obtenir les principaux param?tres g?om?trique de la tail en
%fonction des caract?ristiques des ailes. Cette m?thode est bas?e sur
%l'ouvrage de r?f?rence "Aircraft design, A Systems Engineering Approach"
%de Mohammed H Sadraey et plus pr?cis?ment du chapitre 6.

%Cette m?thode est it?rative et s'appuie sur les ?quilibre statiques
%et dynamique.

%% Data required
c_bar = c_chord;
% - AR
AR = 7; % fixed
% - lambda
lambda = 0.3; % fixed
% - rho
rho = 0.46;
lambda_h = 0.5; % fixed
lambda_v = 0.9; % fixed

Lambda_T = Lambda_LE + 5;

%% Longitudinal trim

% Le coefficient (l/c_bar * S_h/S)=:V_h_bar  est tr?s important et est
% appel? "horizontal tail volume coefficient". Il joue un role tr?s
% important dans la stabilit? longitutinale de l'appareil. Apr?s avoir
% consult? plusieur sources statistiques, ce cefficient sera fix? ? une
% premi?re valeur de 0.6;
V_h_bar = 0.6;
l_arm = 3.3;
% Ainsi, on peut maintenant d?temriner S_h :
S_h = V_h_bar*c_bar*S/l_arm;

K_beta = 0.3*l_cg/l_f + 0.75 *h_f_max/l_f - 0.105;
S_fs = 0.8*l_f*h_f_max;
CN_beta_f = -K_beta*S_fs*l_f/S/b*(1.2)^(1/2); %based on estimation of fuselage geometry
CN_beta_i = 0.012; %because high wing
CN_tot = CN_beta_f + CN_beta_i;
disp(CN_tot);
% voir graphique slide 57
V_v = 0.043; % avec V_v = S_F*l_F/(S*b)
l_F = l_arm; % first guess, distance between cg and fin ac
S_v = V_v*S*b/l_F;

angle = atan(sqrt(S_v/S_h)); % angle de la v_tail
C_L = 0.15*cos(angle);

AR = 3;
S_tail = S_h + S_v;
b_tail = sqrt(AR*S_tail);% span along the tail (one side)
b_h = sin(angle)*b_tail;
b_v = cos(angle)*b_tail/2;
lambda_t = 0.6;
c_root_tail = 2*S_tail/b_tail/(1+lambda_t);
c_tip_tail = lambda_t * c_root_tail;
%Weight
pound = 2.20462262; % kg to lbs
feet = 3.28; % m to ft
[W_mass] = fuel_weight();
W_fuel = W_mass*pound;
W_dg = Mass*pound - 0.45*W_fuel;
N_z = 1.5*3; %Ultimate load factor = 1.5*limit load factor
S = S*feet^2; %[ft^2]
Lambda_W = Lambda_LE*pi/180; % Sweep angle
S_csw = 0.2*S;
q = 8.6292e+03*pound/feet^2;
W_tail_h = 0.016*(N_z*W_dg)^(0.414)*q^(0.168)*(S_tail*feet^2*sin(angle))^(0.896)*...
    (100*0.1/cos(Lambda_W))^(-0.12)*(AR/(cos(Lambda_T))^2)^(0.043)*lambda_t^(-0.02);
Ht = 1;
Hv = 2;
W_tail_v = 0.073*(1+0.2*Ht/Hv)*(N_z*W_dg)^(0.376)*q^(0.122)*(S_tail*feet^2*cos(angle))^0.873...
    *(100*0.1/cos(Lambda_T))^-0.49*(AR/(cos(Lambda_T))^2)^0.357*lambda^0.039;
S = S/feet^2;
fprintf('Tail surface : %.2dft?\n',S_tail);
fprintf('Tail span : %.2dm\n',b_tail);
fprintf('Tail horizontal span : %.2dm\n',b_h);
fprintf('Tail vertical span : %.2dm\n',b_v);
fprintf('Dihedral angle (degrees) : %.2d\n',angle*180/pi);
fprintf('Surface ratio : %.2d\n',S_h/S);
fprintf('Weight of the horizontal tail : %.2dkg\n', W_tail_h/pound);
fprintf('Weight of the vertical tail : %.2dkg\n', W_tail_v/pound);
fprintf('Total weight of the tail : %.2dkg\n', (W_tail_v + W_tail_h)/pound);
W_tail = (W_tail_v + W_tail_h)/pound;