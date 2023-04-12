function [boom, stringers, area] = Boom(cell,stringers, M_x, M_y, M_z, sigma_max,coeff_boom1)

%fonction that computes the min area of the cross section of the booms
%in function of the moments given by the loads of the wing

%voir s'il faut calculer tout ca au tip ou a la root
%et voir si je fais un coeff d'importance relative d'aire ou pas

boom.XZ_up = cell.XZ_up;
boom.XZ_low = cell.XZ_low;

%Centroid computation
centroid.X = (boom.XZ_up(1,1)*coeff_boom1 + boom.XZ_up(1,2) + boom.XZ_low(1,1)*coeff_boom1 + boom.XZ_low(1,2) + sum(stringers.XZ_up1(1,:)) + sum(stringers.XZ_low1(1,:)) + sum(stringers.XZ_up2(1,:)) + sum(stringers.XZ_low2(1,:)))/(stringers.nb+2*coeff_boom1);
centroid.Z = (boom.XZ_up(2,1)*coeff_boom1 + boom.XZ_up(2,2) + boom.XZ_low(2,1)*coeff_boom1 + boom.XZ_low(2,2) + sum(stringers.XZ_up1(2,:)) + sum(stringers.XZ_low1(2,:)) + sum(stringers.XZ_up2(2,:)) + sum(stringers.XZ_low2(2,:)))/(stringers.nb+2*coeff_boom1);

plot(centroid.X ,centroid.Z,'o', 'linewidth', 3,'MarkerSize', 2)

%New system of reference: using the centroid
boom.XZc_up(1,:) = boom.XZ_up(1,:)- centroid.X;
boom.XZc_up(2,:) = boom.XZ_up(2,:)- centroid.Z;
boom.XZc_low(1,:) = boom.XZ_low(1,:)- centroid.X;
boom.XZc_low(2,:) = boom.XZ_low(2,:)- centroid.Z;

stringers.XZc_up1(1,:) = stringers.XZ_up1(1,:)- centroid.X;
stringers.XZc_up1(2,:) = stringers.XZ_up1(2,:)- centroid.Z;

stringers.XZc_low1(1,:) = stringers.XZ_low1(1,:)- centroid.X;
stringers.XZc_low1(2,:) = stringers.XZ_low1(2,:)- centroid.Z;

stringers.XZc_up2(1,:) = stringers.XZ_up2(1,:)- centroid.X;
stringers.XZc_up2(2,:) = stringers.XZ_up2(2,:)- centroid.Z;

stringers.XZc_low2(1,:) = stringers.XZ_low2(1,:)- centroid.X;
stringers.XZc_low2(2,:) = stringers.XZ_low2(2,:)- centroid.Z;


%example: boom.XZc_up(2,1) is the z of the first boom cell, up
%Inertia per unit area
I.xx = boom.XZc_up(2,1)^2*coeff_boom1 + boom.XZc_up(2,2)^2 + boom.XZc_low(2,1)^2*coeff_boom1 + boom.XZc_low(2,2)^2 + sum(stringers.XZc_up1(2,:).^2) + sum(stringers.XZc_low1(2,:).^2) + sum(stringers.XZc_up2(2,:).^2) + sum(stringers.XZc_low2(2,:).^2);
I.xz = boom.XZc_up(2,1)*boom.XZc_up(1,1)*coeff_boom1 + boom.XZc_up(2,2)*boom.XZc_up(1,2) + boom.XZc_low(2,1)*boom.XZc_low(1,1)*coeff_boom1 + boom.XZc_low(2,2)*boom.XZc_low(1,2) + sum(stringers.XZc_up1(2,:).*stringers.XZc_up1(1,:)) + sum(stringers.XZc_low1(2,:).*stringers.XZc_low1(1,:)) + sum(stringers.XZc_up2(2,:).*stringers.XZc_up2(1,:)) + sum(stringers.XZc_low2(2,:).*stringers.XZc_low2(1,:));
I.zz = boom.XZc_up(1,1)^2*coeff_boom1 + boom.XZc_up(1,2)^2 + boom.XZc_low(1,1)^2*coeff_boom1 + boom.XZc_low(1,2)^2 + sum(stringers.XZc_up1(1,:).^2) + sum(stringers.XZc_low1(1,:).^2) + sum(stringers.XZc_up2(1,:).^2) + sum(stringers.XZc_low2(1,:).^2);

%Stress per unit area
boom.sigma_xx_up = sigma_xxx(I,M_x,M_z,boom.XZc_up);
boom.sigma_xx_low = sigma_xxx(I,M_x,M_z,boom.XZc_low);

stringers.sigma_xx_up1 = sigma_xxx(I,M_x,M_z,stringers.XZc_up1);
stringers.sigma_xx_low1 = sigma_xxx(I,M_x,M_z,stringers.XZc_low1);

stringers.sigma_xx_up2 = sigma_xxx(I,M_x,M_z,stringers.XZc_up2);
stringers.sigma_xx_low2 = sigma_xxx(I,M_x,M_z,stringers.XZc_low2);

%Selectioning the worst case
sigma_xx_max = max( [max(abs(boom.sigma_xx_up)) max(abs(boom.sigma_xx_low)) max(abs(stringers.sigma_xx_up1)) max(abs(stringers.sigma_xx_up2)) max(abs(stringers.sigma_xx_low1)) max(abs(stringers.sigma_xx_low2)) ] );
area = sigma_xx_max/sigma_max;

end

function[sigma_xx] = sigma_xxx(I,M_x,M_z,XZc)
    sigma_xx = (-(I.xz*M_x + I.xx*M_z).*XZc(1,:)+(I.zz*M_x + I.xz*M_z).*XZc(2,:))./(I.xx*I.zz - I.xz^2);
end


