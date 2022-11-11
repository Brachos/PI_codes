%Code destin� � obtenir les principaux param�tres g�om�trique de la tail en
%fonction des caract�ristiques des ailes. Cette m�thode est bas�e sur
%l'ouvrage de r�f�rence "Aircraft design, A Systems Engineering Approach"
%de Mohammed H Sadraey et plus pr�cis�ment du chapitre 6.

%Cette m�thode est it�rative et s'appuie sur les �quilibre statiques
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
% Comme pr�sent� sur la figure 6.3 e l'ouvrage et dans l'�quation 6.9, on a
% l'�quation d'�quilibre suivante : M_owf + M_Lwf + M_Lh = 0 
% o� :
% - M_owf est le 'wing/fuselage aerodynamic pitching moment'
% - M_Lwf est le 'lift pitching moment of the wing/fuselage'
% - M_Lh est le 'lift pitching moment of the tail'

% En observant les caract�ristiques g�om�trique de l'avion on peut r��crire
% cette �quation :  M_owf + L_wf*(h*c_bar - h_0*c_bar) - L_h*l_h = 0 (1)
% (voir figure 6.3)

% Apr�s les quelques �tapes expliqu�es dans le livre, on tombe sur
% l'�quation : C_mo_wf + C_L*(h-h_o) - (l/c_bar * S_h/S) * C_Lh = 0 (2)

% Le coefficient (l/c_bar * S_h/S)=:V_h_bar  est tr�s important et est
% appel� "horizontal tail volume coefficient". Il joue un r�le tr�s
% important dans la stabilit� longitutinale de l'appareil. Apr�s avoir
% consult� plusieur sources statistiques, ce cefficient sera fix� � une
% premi�re valeur de 1;
V_h_bar = 1;
% Dans l'�quation (1), on peut trouver les autres coefficients comme suit :
C_mo_wf = C_maf*(AR*cos(gamma)^2)/(AR + 2*cos(gamma)) + 0.01*alpha_t;
% o� :
% - C_maf est le 'wing airfoil section pitching moment coefficient'
% - AR, l'aspect ratio des ailes
% - gamma, le sweep angle des ailes
% - alpha_t, le 'wing twist angle (in degrees)'
C_L = 2*W_avg/(ro*V_c^2*S);
% o� :
% - W_avg est le poids moyen pendant le 'cruising flight'
% - V_c est la vitesse de croisi�re

% L'�quation (2) peut �tre rectifi�e et r��crite en tenant compte de
% l'efficacit� de la tail
% C_mo_wf + C_L*(h-h_o) - eta_h*V_h_bar*C_Lh = 0;
%Nous choisirons une valeur d'efficacit� moyenne de 0.9;
eta_h = 0.9;
% On prend �galement h = 0.25, valeur classique
h = 0.25;
%Admettons que h_o soit connu par l'�tude de l'airfoyl, on peut donc
%trouver C_Lh.
C_Lh = (C_mo_wf + C_L*(h-h_o))/(eta_h*V_h_bar); %(3)
% L'�quation (3) est particuli�rement importante dans le design de la tail
% La longueur l (distance entre le centre a�rodynamique de l'aile et le
% centre a�rdynamique de la tail) est correl� � la taille totale du
% fuselage par le coefficient l/L. Encore une fois, via le tableau 6.2 du
% livre, on donne l/L = 0.5
% Cependant, connaissant le diam�tre maximal du fuselage, S, c_bar et
% V_h_bar, on peut d�terminer la longueur l optimale.
K_c = 1.1; % voir livre page 300)
l = K_c * sqrt(4*c_bar*S*V_h_bar/(pi*D_f));
% Ainsi, on peut maintenant d�temriner S_h :
S_h = V_h_bar*c_bar*S/l;
% L'AR de la tail est d�termin� selon la formule (6.59)
AR_h = 2/3 * AR;