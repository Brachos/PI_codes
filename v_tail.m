function [S_tail, S_h,S_v,c_root_tail, c_tip_tail, angle,l,C_L,Lambda_T, b_tail, b_v, b_h, W_tail] = v_tail(Mass,...
    D_f_max,h_f_max,V_c,c_chord,Lambda_LE,S,l_f,l_cg,b)
%Code destin? ? obtenir les principaux param?tres g?om?trique de la tail en
%fonction des caract?ristiques des ailes. Cette m?thode est bas?e sur
%l'ouvrage de r?f?rence "Aircraft design, A Systems Engineering Approach"
%de Mohammed H Sadraey et plus pr?cis?ment du chapitre 6.

%Cette m?thode est it?rative et s'appuie sur les ?quilibre statiques
%et dynamique.

%% Data required
% - m
% m = 5000;
m = Mass;
% - D_f
% D_f = 0.84;
D_f = D_f_max;
% h_f_max = D_f;
% - V_c
% V_c = 212;
% - alpha_f
% - c_bar
% c_bar = 1.44;
c_bar = c_chord;
% - AR
AR = 7; % fixed
% - lambda
lambda = 0.3; % fixed
% - i_w
% i_w = 0.75;
% - alpha_twist
% alpha_twist = -2;
% - Lambda_LE
% Lambda_LE  = 36.34;
% - Gamma
Gamma = 41; %[?]
% - airfoyl profile
% ac = 0.365;
% - C_L_alpha
% C_L_alpha = 9.8;
% - ro
rho = 0.46;
lambda_h = 0.5; % fixed
lambda_v = 0.9; % ficed
% - C_maf
% - cg_pos
% - h_o
% S = 12.07;
Lambda_T = Lambda_LE + 5;

%% Longitudinal trim
% Comme pr?sent? sur la figure 6.3 e l'ouvrage et dans l'?quation 6.9, on a
% l'?quation d'?quilibre suivante : M_owf + M_Lwf + M_Lh = 0 
% o? :
% - M_owf est le 'wing/fuselage aerodynamic pitching moment'
% - M_Lwf est le 'lift pitching moment of the wing/fuselage'
% - M_Lh est le 'lift pitching moment of the tail'

% En observant les caract?ristiques g?om?trique de l'avion on peut r??crire
% cette ?quation :  M_owf + L_wf*(h*c_bar - h_0*c_bar) - L_h*l_h = 0 (1)
% (voir figure 6.3)

% Apr?s les quelques ?tapes expliqu?es dans le livre, on tombe sur
% l'?quation : C_mo_wf + C_L*(h-h_o) - (l/c_bar * S_h/S) * C_Lh = 0 (2)


% Dans l'?quation (1), on peut trouver les autres coefficients comme suit :
% C_mo_wf = C_maf*(AR*cos(Alpha)^2)/(AR + 2*cos(Alpha)) + 0.01*alpha_t;
% o? :
% - C_maf est le 'wing airfoil section pitching moment coefficient'
% - AR, l'aspect ratio des ailes
% - gamma, le sweep angle des ailes
% - alpha_t, le 'wing twist angle (in degrees)'
% C_L = 2*W_avg/(ro*V_c^2*S);
% o? :
% - W_avg est le poids moyen pendant le 'cruising flight'
% - V_c est la vitesse de croisi?re

% L'?quation (2) peut ?tre rectifi?e et r??crite en tenant compte de
% l'efficacit? de la tail
% C_mo_wf + C_L*(h-h_o) - eta_h*V_h_bar*C_Lh = 0;
%Nous choisirons une valeur d'efficacit? moyenne de 0.9;
eta_h = 0.9;
% On prend ?galement h = 0.25, valeur classique
% h = 0.25;

% Le coefficient (l/c_bar * S_h/S)=:V_h_bar  est tr?s important et est
% appel? "horizontal tail volume coefficient". Il joue un role tr?s
% important dans la stabilit? longitutinale de l'appareil. Apr?s avoir
% consult? plusieur sources statistiques, ce cefficient sera fix? ? une
% premi?re valeur de 0.75;
V_h_bar = 0.6;
% K_c = 1.1; % voir livre page 300)
% l_opt = K_c * sqrt(4*c_bar*S*V_h_bar/(pi*D_f));
% disp(l_opt);
% Ratio des longueur (voir page 276)
l_ratio = 0.45; %0.55 un peu trop ?lev?
l = l_f*l_ratio;
l = l_f - l_cg - 0.5;
l = 3.3;
% Ainsi, on peut maintenant d?temriner S_h :
S_h = V_h_bar*c_bar*S/l;

% fprintf('Lift coefficient of the tail = %d\n',C_L);

% Selon le profil de l'aile, on a le pourcentage de la longueur auquel se
% trouve le aerodynamic center (ex. 32% of the MAC)
% Selon le centre de gravit? du fuselage, en pourcentage de la longueur du
% fuselage, on peut trouver X_apex
% cg_pos est la poition du centre de gravit? par rapport ? la position au
% centre a?rodynamique
% cg_pos = 2;
%X_apex = -h_o*c_bar + 0.32*L_f + cg_pos; %!!!!!!! revoir les valeurs 0.23 et 0.32 !!!!
%X_cg = h_o*c_bar - cg_pos;
%h = X_cg/c_bar;

%Admettons que h_o soit connu par l'?tude de l'airfoyl, on peut donc
%trouver C_Lh.
%C_Lh = (C_mo_wf + C_L*(h-h_o))?/(eta_h*V_h_bar); %(3)
% L'?quation (3) est particuli?rement importante dans le design de la tail
% La longueur l (distance entre le centre a?rodynamique de l'aile et le
% centre a?rdynamique de la tail) est correl? ? la taille totale du
% fuselage par le coefficient l/L. Encore une fois, via le tableau 6.2 du
% livre, on donne l/L = 0.5
% Cependant, connaissant le diam?tre maximal du fuselage, S, c_bar et
% V_h_bar, on peut d?terminer la longueur l optimale.

% L'AR de la tail est d?termin? selon la formule (6.59)
AR_h = 2/3 * AR;
b_h = sqrt(S_h*AR_h);

K_beta = 0.3*l_cg/l_f + 0.75 *h_f_max/l_f - 0.105;
S_fs = 0.8*l_f*h_f_max;
CN_beta_f = -K_beta*S_fs*l_f/S/b;
CN_beta_i = 0.012; %because mid wing
CN_tot = CN_beta_f + CN_beta_i;
disp(CN_tot);
% voir graphique slide 57
V_v = 0.04; % avec V_v = S_F*l_F/(S*b)
l_F = l; % first guess, distance between cg and fin ac
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
Lambda_W = 36.34*pi/180; % Sweep angle
S_csw = 0.2*S;
q = 8.6292e+03*pound/feet^2;
W_tail_h = 0.016*(N_z*W_dg)^(0.414)*q^(0.168)*(S_tail*feet^2*sin(angle))^(0.896)*...
    (100*0.1/cos(Lambda_W))^(-0.12)*(AR/(cos(Lambda_T))^2)^(0.043)*lambda_t^(-0.02);
%Ht/Hv ??? = 0.5 because no information about V-tail
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