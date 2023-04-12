function [] = skin_size(boom_root,boom_tip, stringers_root, stringers_tip, M_x, M_z)
%function to compute the size of the skin

%inertia with respect to the centroid
I.xx = boom.XZc_up(2,1)^2*boom.Area(1) + boom.XZc_up(2,2)^2*boom.Area(2) + boom.XZc_low(2,1)^2*boom.Area(1) + boom.XZc_low(2,2)^2*boom.Area(2) + sum(stringers_root.XZc_up1(2,:).^2*stringers_root.Area) + sum(stringers_root.XZc_low1(2,:).^2*stringers_root.Area) + sum(stringers_root.XZc_up2(2,:).^2*stringers_root.Area) + sum(stringers_root.XZc_low2(2,:).^2*stringers_root.Area);
I.xz = boom.XZc_up(2,1)*boom.XZc_up(1,1)*boom.Area(1) + boom.XZc_up(2,2)*boom.XZc_up(1,2)*boom.Area(2) + boom.XZc_low(2,1)*boom.XZc_low(1,1)*coeff_boom1*boom.Area(1) + boom.XZc_low(2,2)*boom.XZc_low(1,2)*boom.Area(2) + sum(stringers_root.XZc_up1(2,:).*stringers_root.XZc_up1(1,:)*stringers_root.Area) + sum(stringers_root.XZc_low1(2,:).*stringers_root.XZc_low1(1,:)*stringers_root.Area) + sum(stringers_root.XZc_up2(2,:).*stringers_root.XZc_up2(1,:)*stringers_root.Area) + sum(stringers_root.XZc_low2(2,:).*stringers_root.XZc_low2(1,:)*stringers_root.Area);
I.zz = boom.XZc_up(1,1)^2*boom.Area(1) + boom.XZc_up(1,2)^2*boom.Area(2) + boom.XZc_low(1,1)^2*boom.Area(1) + boom.XZc_low(1,2)^2*boom.Area(2) + sum(stringers_root.XZc_up1(1,:).^2*stringers_root.Area) + sum(stringers_root.XZc_low1(1,:).^2*stringers_root.Area) + sum(stringers_root.XZc_up2(1,:).^2*stringers_root.Area) + sum(stringers_root.XZc_low2(1,:).^2*stringers_root.Area);

%Direct stress computation of the booms and of the stringers
boom_root.sigma_xx_up = (-(I.xz*M_x + I.xx*M_z).*boom_root.XZc_up(1,:)+(I.zz*M_x + I.xz*M_z).*boom_root.XZc_up(2,:))./(I.xx*I.zz - I.xz^2);
boom_root.sigma_xx_low =  (-(I.xz*M_x + I.xx*M_z).*boom_root.XZc_low(1,:)+(I.zz*M_x + I.xz*M_z).*boom_root.XZc_low(2,:))./(I.xx*I.zz - I.xz^2);


stringers_root.sigma_xx_up1 = (-(I.xz*M_x + I.xx*M_z).*stringers.XZc_up1(1,:)+(I.zz*M_x + I.xz*M_z).*stringers.XZc_up1(2,:))./(I.xx*I.zz - I.xz^2);
stringers_root.sigma_xx_up2 = (-(I.xz*M_x + I.xx*M_z).*stringers.XZc_up2(1,:)+(I.zz*M_x + I.xz*M_z).*stringers.XZc_up2(2,:))./(I.xx*I.zz - I.xz^2);
   
stringers_root.sigma_xx_low1 = (-(I.xz*M_x + I.xx*M_z).*stringers.XZc_low1(1,:)+(I.zz*M_x + I.xz*M_z).*stringers.XZc_low1(2,:))./(I.xx*I.zz - I.xz^2);
stringers_root.sigma_xx_low2 = (-(I.xz*M_x + I.xx*M_z).*stringers.XZc_low2(1,:)+(I.zz*M_x + I.xz*M_z).*stringers.XZc_low2(2,:))./(I.xx*I.zz - I.xz^2);

end

