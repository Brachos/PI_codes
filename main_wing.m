function [] = main_wing()
%MAIN_WING 
%input for all the functions computing the structure of the wing
[cell_root, cell_tip, stringers_root, stringers_tip] = wing_geom();

sigma_max = ....

%a voir les arguments qu'on lui donne, ou comment on le lance!!
%pcq j'aimerais que en sortie j'ai des vecteurs Mx,My et Mz
%soit on modifie les codes et on renvoie le vecteur
%soit on forme le vecteur dans ce main ci
%soit on lance avant la fonction loads, qui va renvoyer le bon vecteur!!

[Tx,Ty,Tz,Mx,My,Mz] = wing_load(W,Wing_loading,i,y_cg,y_ac,Mom_wing,n,Wing);


%--------------------------------
%Boom area computation

%il faut voir si on garde ca ou si on ajoute aussi des coeff pour les
%autres stringers et tout
coeff_boom1 = 10; %to have an area of the booms of the 1st cell 10 times bigger


%voir s'il faut calculer tout ca au tip ou a la root
%et voir si je fais un coeff d'importance relative d'aire ou pas
for m = 1 : length(Mx)
    %tous les endroits de la maneuver
[boom_root, boom_tip stringers_root, area_min(m)] = Boom(cell_root,cell_tip,stringers_root, Mx(m), My(m), Mz(m), sigma_max, coeff_boom1);
end

Area = max(area_min);

boom_root.Area(1) = Area*coeff_boom1;
boom_root.Area(2) = Area;
stringers_root.Area = Area;


%------------------------------------------------
%skin



end

