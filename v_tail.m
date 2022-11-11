%Code destiné à obtenir les principaux paramètres géométrique de la tail en
%fonction des caractéristiques des ailes. Cette méthode est basée sur
%l'ouvrage de référence "Aircraft design, A Systems Engineering Approach"
%de Mohammed H Sadraey et plus précisément du chapitre 6.

%Cette méthode est itérative et s'appuie sur les équilibre statiques
%et dynamique.

%% Data required
% - c_bar
% - AR
% - S
% - ro
% - h_o ~ [0.2;0.25]
% - C_maf
% - gamma
% - alpha_t
% - W_avg
% - V_c
% - D_f

%% Longitudinal trim
% Comme présenté sur la figure 6.3 e l'ouvrage et dans l'équation 6.9, on a
% l'équation d'équilibre suivante : M_owf + M_Lwf + M_Lh = 0 
% où :
% - M_owf est le 'wing/fuselage aerodynamic pitching moment'
% - M_Lwf est le 'lift pitching moment of the wing/fuselage'
% - M_Lh est le 'lift pitching moment of the tail'

% En observant les caractéristiques géométrique de l'avion on peut réécrire
% cette équation :  M_owf + L_wf*(h*c_bar - h_0*c_bar) - L_h*l_h = 0 (1)
% (voir figure 6.3)

% Après les quelques étapes expliquées dans le livre, on tombe sur
% l'équation : C_mo_wf + C_L*(h-h_o) - (l/c_bar * S_h/S) * C_Lh = 0 (2)

% Le coefficient (l/c_bar * S_h/S)=:V_h_bar  est très important et est
% appelé "horizontal tail volume coefficient". Il joue un rôle très
% important dans la stabilité longitutinale de l'appareil. Après avoir
% consulté plusieur sources statistiques, ce cefficient sera fixé à une
% première valeur de 1;
V_h_bar = 1;
% Dans l'équation (1), on peut trouver les autres coefficients comme suit :
C_mo_wf = C_maf*(AR*cos(gamma)^2)/(AR + 2*cos(gamma)) + 0.01*alpha_t;
% où :
% - C_maf est le 'wing airfoil section pitching moment coefficient'
% - AR, l'aspect ratio des ailes
% - gamma, le sweep angle des ailes
% - alpha_t, le 'wing twist angle (in degrees)'
C_L = 2*W_avg/(ro*V_c^2*S);
% où :
% - W_avg est le poids moyen pendant le 'cruising flight'
% - V_c est la vitesse de croisière

% L'équation (2) peut être rectifiée et réécrite en tenant compte de
% l'efficacité de la tail
% C_mo_wf + C_L*(h-h_o) - eta_h*V_h_bar*C_Lh = 0;
%Nous choisirons une valeur d'efficacité moyenne de 0.9;
eta_h = 0.9;
% On prend également h = 0.25, valeur classique
h = 0.25;
%Admettons que h_o soit connu par l'étude de l'airfoyl, on peut donc
%trouver C_Lh.
C_Lh = (C_mo_wf + C_L*(h-h_o))/(eta_h*V_h_bar); %(3)
% L'équation (3) est particulièrement importante dans le design de la tail
% La longueur l (distance entre le centre aérodynamique de l'aile et le
% centre aérdynamique de la tail) est correlé à la taille totale du
% fuselage par le coefficient l/L. Encore une fois, via le tableau 6.2 du
% livre, on donne l/L = 0.5
% Cependant, connaissant le diamètre maximal du fuselage, S, c_bar et
% V_h_bar, on peut déterminer la longueur l optimale.
K_c = 1.1; % voir livre page 300)
l = K_c * sqrt(4*c_bar*S*V_h_bar/(pi*D_f));
% Ainsi, on peut maintenant détemriner S_h :
S_h = V_h_bar*c_bar*S/l;
% L'AR de la tail est déterminé selon la formule (6.59)
AR_h = 2/3 * AR;