    function [cell_root, cell_tip, stringers_root, stringers_tip] = wing_geom()
    %WING_GEOM describes the geometry of the wing including the stringers and
    %the booms
    %à voir plus tard:
    %-checker mon code pour la tipe, voir avec le twist + la translation!
    %-est ce que je changerais pas les cell en booms puis après poru faire les
    %cels je les relie juste?
    %sinon ca fait vrmt comme le code de l'année passée

    %After that compute the reste to have the centroid
    %after the boom!!!

    [b,S,CL_alpha,CD_alpha,CL,CD,D,c_root,c_tip,c_AC,x_AC,y_AC,V_fuel,sweep,c,alpha_L0,tap,cl_alpha,theta_tip, AR] = wing(0.7,30000, 4236,2.5)
    close all
    coord   = dlmread('SC(2)-0714.txt');
    xz_up = coord(1:103,:);
    xz_low = coord(104:end,:);

    %% At the root
    airf_root.XZ_up = c_root*xz_up; %position along x and z of the airfoil's upper part at the root
    airf_root.XZ_low = c_root*xz_low; %lower part of the airfoil

    plot(airf_root.XZ_up(:,1), airf_root.XZ_up(:,2),'linewidth', 2, 'MarkerSize', 11')
    hold on
    plot(airf_root.XZ_low(:,1),airf_root.XZ_low(:,2),'linewidth', 2, 'MarkerSize', 11')
    axis equal

    %first cell (on the left)
    %position: at the center of gravity of the airfoil's root
    cell_root.x(1) = c_root/4;
    [cell_root] = Cell_coord(airf_root, cell_root, 1);

    %second cell (on the right)
    %position: +/- at the begining of the flaps
    cell_root.x(2) = 0.7*c_root;
    [cell_root] = Cell_coord(airf_root, cell_root, 2);

    plot([cell_root.XZ_up(1,1) cell_root.XZ_low(1,1)], [cell_root.XZ_up(2,1) cell_root.XZ_low(2,1)],'linewidth', 2)
    plot([cell_root.XZ_up(1,2) cell_root.XZ_low(1,2)], [cell_root.XZ_up(2,2) cell_root.XZ_low(2,2)],'linewidth', 2)

    plot([cell_root.XZ_up(1,1) cell_root.XZ_low(1,1)], [cell_root.XZ_up(2,1) cell_root.XZ_low(2,1)],'o', 'linewidth', 5,'MarkerSize', 2)
    plot([cell_root.XZ_up(1,2) cell_root.XZ_low(1,2)], [cell_root.XZ_up(2,2) cell_root.XZ_low(2,2)],'o', 'linewidth', 5,'MarkerSize', 2)

    %stringers on the first cell
    nb_str_root_1 = 4;
    [stringers_root.XZ_up1] = stringers_coord(airf_root.XZ_up, nb_str_root_1, 1, cell_root.XZ_up);
    stringers_0 = [0; 0];
    stringers_root.XZ_up1 = [stringers_0  stringers_root.XZ_up1];
    [stringers_root.XZ_low1] = stringers_coord(airf_root.XZ_low, nb_str_root_1, 1, cell_root.XZ_low);

    %stringers on the second cell
    nb_str_root_2 = 7;
    [stringers_root.XZ_up2] = stringers_coord(airf_root.XZ_up, nb_str_root_2, 2,  cell_root.XZ_up);
    [stringers_root.XZ_low2] = stringers_coord(airf_root.XZ_low, nb_str_root_2, 2, cell_root.XZ_low);

    plot(stringers_root.XZ_up1(1,:), stringers_root.XZ_up1(2,:),'o', 'linewidth', 3,'MarkerSize', 2)
    plot(stringers_root.XZ_low1(1,:),stringers_root.XZ_low1(2,:),'o', 'linewidth', 3,'MarkerSize', 2)

    plot(stringers_root.XZ_up2(1,:), stringers_root.XZ_up2(2,:),'o', 'linewidth', 3,'MarkerSize', 2)
    plot(stringers_root.XZ_low2(1,:),stringers_root.XZ_low2(2,:),'o', 'linewidth', 3,'MarkerSize', 2)

    stringers_root.nb = nb_str_root_2 * 2 + nb_str_root_1 * 2 +5;


    %% At the tip

    %have to take into account the twist angle!! /!\
    %AND the translation!!!
    airf_tip.XZ_up = c_tip*xz_up; %position along x and Z of the airfoil's upper part at the tip
    airf_tip.XZ_low = c_tip*xz_low; %lower part of the airfoil


    plot(airf_tip.XZ_up(:,1), airf_tip.XZ_up(:,2),'linewidth', 2)
    hold on
    plot(airf_tip.XZ_low(:,1),airf_tip.XZ_low(:,2),'linewidth', 2)
    axis equal


    %third cell (on the left)
    %position: at the center of gravity of the airfoil's tip
    cell_tip.x(1) = c_tip/4;
    [cell_tip] = Cell_coord(airf_tip, cell_tip, 1);

    %fourth cell (on the right)
    %position: +/- at the begining of the flaps
    cell_tip.x(2) = 0.7*c_tip;
    [cell_tip] = Cell_coord(airf_tip, cell_tip, 2);

    plot([cell_tip.XZ_up(1,1) cell_tip.XZ_low(1,1)], [cell_tip.XZ_up(2,1) cell_tip.XZ_low(2,1)],'linewidth', 2)
    plot([cell_tip.XZ_up(1,2) cell_tip.XZ_low(1,2)], [cell_tip.XZ_up(2,2) cell_tip.XZ_low(2,2)],'linewidth', 2)

    plot([cell_tip.XZ_up(1,1) cell_tip.XZ_low(1,1)], [cell_tip.XZ_up(2,1) cell_tip.XZ_low(2,1)],'o', 'linewidth', 5,'MarkerSize', 2)
    plot([cell_tip.XZ_up(1,2) cell_tip.XZ_low(1,2)], [cell_tip.XZ_up(2,2) cell_tip.XZ_low(2,2)],'o', 'linewidth', 5,'MarkerSize', 2)

    %stringers on the first cell
    nb_str_tip_1 = nb_str_root_1;
    [stringers_tip.XZ_up1] = stringers_coord(airf_tip.XZ_up, nb_str_tip_1, 1, cell_tip.XZ_up);
    stringers_0_tip = airf_tip.XZ_up(1,:);
    stringers_tip.XZ_up1 = [stringers_0_tip'  stringers_tip.XZ_up1];
    [stringers_tip.XZ_low1] = stringers_coord(airf_tip.XZ_low, nb_str_tip_1, 1, cell_tip.XZ_low);

    %stringers on the second cell
    nb_str_tip_2 = nb_str_root_2;
    [stringers_tip.XZ_up2] = stringers_coord(airf_tip.XZ_up, nb_str_tip_2, 2, cell_tip.XZ_up);
    [stringers_tip.XZ_low2] = stringers_coord(airf_tip.XZ_low, nb_str_tip_2, 2, cell_tip.XZ_low);
    %[w, j, l] = Boom(cell_root,stringers_root,1,1,1,10,10)

    plot(stringers_tip.XZ_up1(1,:), stringers_tip.XZ_up1(2,:),'o', 'linewidth', 3,'MarkerSize', 2)
    plot(stringers_tip.XZ_low1(1,:),stringers_tip.XZ_low1(2,:),'o', 'linewidth', 3,'MarkerSize', 2)

    plot(stringers_tip.XZ_up2(1,:), stringers_tip.XZ_up2(2,:),'o', 'linewidth', 3,'MarkerSize', 2)
    plot(stringers_tip.XZ_low2(1,:),stringers_tip.XZ_low2(2,:),'o', 'linewidth', 3,'MarkerSize', 2)

    end
    %%

    function [cell] = Cell_coord(airf, cell, num_cell)

    [~, cell.index_up(num_cell)] = min(abs(airf.XZ_up(:,1)-cell.x(num_cell)));
    cell.XZ_up(1,num_cell) = airf.XZ_up(cell.index_up(num_cell),1);
    cell.XZ_up(2,num_cell) = airf.XZ_up(cell.index_up(num_cell),2);

    [~, cell.index_low(num_cell)] = min(abs(airf.XZ_low(:,1)-cell.x(num_cell)));
    cell.XZ_low(1,num_cell) = airf.XZ_low(cell.index_low(num_cell),1);
    cell.XZ_low(2,num_cell) = airf.XZ_low(cell.index_low(num_cell),2);

    end

    %la fonction doit etre appellée 4x si mes calculs sont bons!!
    %donc deux fois pour le up, 2 fois pour le low
    %et puis encore x2 si on prend le tip et le root

    function[XZ] = stringers_coord(airf, nb_str, num_cell, cell)
    %airf = the x and Z coord of the airfoil
    if num_cell == 1
        dx = (cell(1,1)-airf(1,1))/(nb_str+1);
        x_str = ones(1,nb_str);
        for n = 1 : nb_str
            x_str(n) = airf(1,1) + dx*n;
        end
    end
    if num_cell == 2
        dx = (cell(1,2)-cell(1,1))/(nb_str+1);
        x_str = ones(1,nb_str);
        for n = 1 : nb_str
            x_str(n) = cell(1,1) + dx*n;
        end
    end
    %on prend en argument la courbe, et la coordonnée sur l'axe x, il faut
    %enregistrer les coord x_Z et renvoyer ca
    %index =
    for n = 1 : nb_str
        [~, index(n)] = min(abs(airf(:,1)-x_str(n)));
        XZ(1,n) = airf(index(n),1);
        XZ(2,n) = airf(index(n),2);
    end
    end



