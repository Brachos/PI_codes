function [thickness] = skin_size(boom_root,boom_tip, stringers_root, stringers_tip, M_x, M_z, T_x, T_z,airf_root,mu,mu_ref,tau_max)
%function to compute the size of the skin

%inertia with respect to the centroid
%/!\ est ce que c'est boom_tip ou boom_root?
I.xx = boom_root.XZc_up(2,1)^2*boom_root.Area(1) + boom_root.XZc_up(2,2)^2*boom_root.Area(2) + boom_root.XZc_low(2,1)^2*boom_root.Area(1) + boom_root.XZc_low(2,2)^2*boom_root.Area(2) + sum(stringers_root.XZc_up1(2,:).^2*stringers_root.Area) + sum(stringers_root.XZc_low1(2,:).^2*stringers_root.Area) + sum(stringers_root.XZc_up2(2,:).^2*stringers_root.Area) + sum(stringers_root.XZc_low2(2,:).^2*stringers_root.Area);
I.xz = boom_root.XZc_up(2,1)*boom_root.XZc_up(1,1)*boom_root.Area(1) + boom_root.XZc_up(2,2)*boom_root.XZc_up(1,2)*boom_root.Area(2) + boom_root.XZc_low(2,1)*boom_root.XZc_low(1,1)*boom_root.Area(1) + boom_root.XZc_low(2,2)*boom_root.XZc_low(1,2)*boom_root.Area(2) + sum(stringers_root.XZc_up1(2,:).*stringers_root.XZc_up1(1,:)*stringers_root.Area) + sum(stringers_root.XZc_low1(2,:).*stringers_root.XZc_low1(1,:)*stringers_root.Area) + sum(stringers_root.XZc_up2(2,:).*stringers_root.XZc_up2(1,:)*stringers_root.Area) + sum(stringers_root.XZc_low2(2,:).*stringers_root.XZc_low2(1,:)*stringers_root.Area);
I.zz = boom_root.XZc_up(1,1)^2*boom_root.Area(1) + boom_root.XZc_up(1,2)^2*boom_root.Area(2) + boom_root.XZc_low(1,1)^2*boom_root.Area(1) + boom_root.XZc_low(1,2)^2*boom_root.Area(2) + sum(stringers_root.XZc_up1(1,:).^2*stringers_root.Area) + sum(stringers_root.XZc_low1(1,:).^2*stringers_root.Area) + sum(stringers_root.XZc_up2(1,:).^2*stringers_root.Area) + sum(stringers_root.XZc_low2(1,:).^2*stringers_root.Area);

%Direct stress computation of the booms and of the stringers
%the direct stress in along y in the syst of reference of the airplane
%(along the wing span)
boom_root.sigma_xx_up = (-(I.xz*M_x + I.xx*M_z).*boom_root.XZc_up(1,:)+(I.zz*M_x + I.xz*M_z).*boom_root.XZc_up(2,:))./(I.xx*I.zz - I.xz^2);
boom_root.sigma_xx_low =  (-(I.xz*M_x + I.xx*M_z).*boom_root.XZc_low(1,:)+(I.zz*M_x + I.xz*M_z).*boom_root.XZc_low(2,:))./(I.xx*I.zz - I.xz^2);


stringers_root.sigma_xx_up1 = (-(I.xz*M_x + I.xx*M_z).*stringers_root.XZc_up1(1,:)+(I.zz*M_x + I.xz*M_z).*stringers_root.XZc_up1(2,:))./(I.xx*I.zz - I.xz^2);
stringers_root.sigma_xx_up2 = (-(I.xz*M_x + I.xx*M_z).*stringers_root.XZc_up2(1,:)+(I.zz*M_x + I.xz*M_z).*stringers_root.XZc_up2(2,:))./(I.xx*I.zz - I.xz^2);
   
stringers_root.sigma_xx_low1 = (-(I.xz*M_x + I.xx*M_z).*stringers_root.XZc_low1(1,:)+(I.zz*M_x + I.xz*M_z).*stringers_root.XZc_low1(2,:))./(I.xx*I.zz - I.xz^2);
stringers_root.sigma_xx_low2 = (-(I.xz*M_x + I.xx*M_z).*stringers_root.XZc_low2(1,:)+(I.zz*M_x + I.xz*M_z).*stringers_root.XZc_low2(2,:))./(I.xx*I.zz - I.xz^2);

%% Delta x and z due to tapering

boom.delta_x_up = (boom_tip.XZ_up(1,:)' - boom_root.XZ_up(1,:)')./(boom_tip.Y-boom_root.Y);
boom.delta_x_low = (boom_tip.XZ_low(1,:)' - boom_root.XZ_low(1,:)')./(boom_tip.Y-boom_root.Y);

boom.delta_z_up = (boom_tip.XZ_up(2,:)' - boom_root.XZ_up(2,:)')./(boom_tip.Y-boom_root.Y);
boom.delta_z_low = (boom_tip.XZ_low(2,:)' - boom_root.XZ_low(2,:)')./(boom_tip.Y-boom_root.Y);


stringers.delta_x_up1 = (stringers_tip.XZ_up1(1,:) - stringers_root.XZ_up1(1,:))/(boom_tip.Y(1)-boom_root.Y(1));
stringers.delta_x_up2 =(stringers_tip.XZ_up2(1,:) - stringers_root.XZ_up2(1,:))/(boom_tip.Y(1)-boom_root.Y(1));
stringers.delta_x_low1 = (stringers_tip.XZ_low1(1,:) - stringers_root.XZ_low1(1,:))/(boom_tip.Y(1)-boom_root.Y(1));
stringers.delta_x_low2 = (stringers_tip.XZ_low2(1,:) - stringers_root.XZ_low2(1,:))/(boom_tip.Y(1)-boom_root.Y(1));

stringers.delta_z_up1 = (stringers_tip.XZ_up1(2,:) - stringers_root.XZ_up1(2,:))/(boom_tip.Y(1)-boom_root.Y(1));
stringers.delta_z_up2 = (stringers_tip.XZ_up2(2,:) - stringers_root.XZ_up2(2,:))/(boom_tip.Y(1)-boom_root.Y(1));
stringers.delta_z_low1 = (stringers_tip.XZ_low1(2,:) - stringers_root.XZ_low1(2,:))/(boom_tip.Y(1)-boom_root.Y(1));
stringers.delta_z_low2 = (stringers_tip.XZ_low2(2,:) - stringers_root.XZ_low2(2,:))/(boom_tip.Y(1)-boom_root.Y(1));
 
%% Stress in the web

boom.P_y_up = boom_root.sigma_xx_up.*boom_root.Area;
boom.P_x_up = boom.P_y_up.*boom.delta_x_up';
boom.P_z_up = boom.P_y_up.*boom.delta_z_up';

boom.P_y_low = boom_root.sigma_xx_low.*boom_root.Area;
boom.P_x_low = boom.P_y_low.*boom.delta_x_low';
boom.P_z_low = boom.P_y_low.*boom.delta_z_low';

stringers.P_y_up1 = stringers_root.sigma_xx_up1*stringers_root.Area;
stringers.P_x_up1 = stringers.P_y_up1.*stringers.delta_x_up1;
stringers.P_z_up1 = stringers.P_y_up1.*stringers.delta_z_up1;

stringers.P_y_up2 = stringers_root.sigma_xx_up2*stringers_root.Area;
stringers.P_x_up2 = stringers.P_y_up2.*stringers.delta_x_up2;
stringers.P_z_up2 = stringers.P_y_up2.*stringers.delta_z_up2;

stringers.P_y_low1 = stringers_root.sigma_xx_low1*stringers_root.Area;
stringers.P_x_low1 = stringers.P_y_low1.*stringers.delta_x_low1;
stringers.P_z_low1 = stringers.P_y_low1.*stringers.delta_z_low1;

stringers.P_y_low2 = stringers_root.sigma_xx_low2*stringers_root.Area;
stringers.P_x_low2 = stringers.P_y_low2.*stringers.delta_x_low2;
stringers.P_z_low2 = stringers.P_y_low2.*stringers.delta_z_low2;


T_webb.x = T_x -(sum(boom.P_x_up) + sum(boom.P_x_low) + sum(stringers.P_x_up1) + sum(stringers.P_x_up2) + sum(stringers.P_x_low1) + sum(stringers.P_x_low2));
T_webb.z = T_z -(sum(boom.P_z_up) + sum(boom.P_z_low) + sum(stringers.P_z_up1) + sum(stringers.P_z_up2) + sum(stringers.P_z_low1) + sum(stringers.P_z_low2));


%% Open shear flux 

%first segment: from 4 to 1: curved 
nb_str = length(stringers_root.XZ_up2(1,:));
q0_41c = zeros(1,nb_str+1); 
q0_41c(1) = 0;

for j=1:nb_str
    q0_41c(j+1) = q0_41c(j) - (I.zz*T_webb.z - I.xz*T_webb.x) / (I.xx*I.zz - I.xz^2) * sum(stringers_root.XZc_up2(2,nb_str-j+1) * stringers_root.Area);
    q0_41c(j+1) = q0_41c(j+1) -(I.xx *T_webb.x - I.xz*T_webb.z)/(I.xx*I.zz - I.xz^2) * sum(stringers_root.XZc_up2(1,nb_str-j+1) * stringers_root.Area);
end


%second segment: from 1 to 2: curved
%the parcours is done anticlockwise, wich is the oposite order of the upper
%coordinates and the right order of the lower coordinates
nb_str1 =  length(stringers_root.XZ_up1(1,:));
nb_str2 = length(stringers_root.XZ_low1(1,:));
nb_str = nb_str1 + nb_str2;
q0_12c = zeros(1,nb_str+1);
q0_12c(1) = 0; %q0_41c(end);
%upper coordinates
for j=1:nb_str1
    q0_12c(j+1) = q0_12c(j) - (I.zz*T_webb.z - I.xz*T_webb.x) / (I.xx*I.zz - I.xz^2) * (stringers_root.XZc_up1(2,nb_str1-j+1) * stringers_root.Area);
    q0_12c(j+1) = q0_12c(j+1) -(I.xx *T_webb.x - I.xz*T_webb.z)/(I.xx*I.zz - I.xz^2) * (stringers_root.XZc_up1(1,nb_str1-j+1) * stringers_root.Area);
end
%lower coordinates
for j=1:nb_str2
    q0_12c(j+1+nb_str1) = q0_12c(j+nb_str1) - (I.zz*T_webb.z - I.xz*T_webb.x) / (I.xx*I.zz - I.xz^2) * (stringers_root.XZc_low1(2,j) * stringers_root.Area);
    q0_12c(j+1+nb_str1) = q0_12c(j+1+nb_str1) -(I.xx *T_webb.x - I.xz*T_webb.z)/(I.xx*I.zz - I.xz^2) * (stringers_root.XZc_low1(1,j) * stringers_root.Area);
end


%third segment: from 1 to 2: straight line
q0_12s = q0_41c(end) - (I.zz*T_webb.z - I.xz*T_webb.x) / (I.xx*I.zz - I.xz^2) * (boom_root.XZc_up(2,1) * boom_root.Area(1));
q0_12s = q0_12s - (I.xx *T_webb.x - I.xz*T_webb.z)/(I.xx*I.zz - I.xz^2) * (boom_root.XZc_up(1,1) * boom_root.Area(1));


%fourth segment: from 2 to 3: curve
nb_str = length(stringers_root.XZ_low2(1,:));
q0_23c = zeros(1,nb_str+1); 
q0_23c(1) = q0_12s + q0_12c(end) - (I.zz*T_webb.z - I.xz*T_webb.x) / (I.xx*I.zz - I.xz^2) * (boom_root.XZc_low(2,1) * boom_root.Area(1));
q0_23c(1) = q0_23c(1) - (I.xx *T_webb.x - I.xz*T_webb.z)/(I.xx*I.zz - I.xz^2) * (boom_root.XZc_low(1,1) * boom_root.Area(1));

for j=1:nb_str
    q0_23c(j+1) = q0_23c(j) - (I.zz*T_webb.z - I.xz*T_webb.x) / (I.xx*I.zz - I.xz^2) * sum(stringers_root.XZc_low2(2,j) * stringers_root.Area);
    q0_23c(j+1) = q0_23c(j+1) -(I.xx *T_webb.x - I.xz*T_webb.z)/(I.xx*I.zz - I.xz^2) * sum(stringers_root.XZc_low2(1,j) * stringers_root.Area);
end


%fifth segment: from 3 to 4: straight
q0_34s = q0_23c(end) - (I.zz*T_webb.z - I.xz*T_webb.x) / (I.xx*I.zz - I.xz^2) * (boom_root.XZc_low(2,2) * boom_root.Area(2));
q0_34s = q0_34s - (I.xx *T_webb.x - I.xz*T_webb.z)/(I.xx*I.zz - I.xz^2) * (boom_root.XZc_up(1,2) * boom_root.Area(2));


%% adimensional length

%first segment: from 1 to 4: curved 

%/!\voir s'il faut enregistrer d'autres données!! ou faire une structure

L41c_index = [boom_root.index_up(1) stringers_root.index_up2 boom_root.index_up(2)];

for i = 1:length(L41c_index)-1
    index_seg = [L41c_index(i) L41c_index(i+1)];
    x = airf_root.XZ_up(index_seg(1) : index_seg(2), 1);
    y = airf_root.XZ_up(index_seg(1) : index_seg(2), 2);
    L=0; 
    A=0;
    for j = 1: length(x)-1
        L =  L + sqrt((x(j) - x(j+1))^2 + (y(j) - y(j+1))^2);
        
        h(j) = abs( (x(j+1) - x(j))*y(j) - x(j)*(y(j+1) - y(j)) )/( (x(j+1) - x(j))^2 + (y(j+1) - y(j))^2 )^(1/2);
        A = A + h(j)*1/2 * sqrt((x(j) - x(j+1))^2 + (y(j) - y(j+1))^2);
    end
    L41c_l(i) = L;
    L41c_A(i) = A;
end

L41c_l = flip(L41c_l);
L41c_A = flip(L41c_A);
L41c_L = sum(L41c_l);
L41c_lbar_t = L41c_l * mu_ref/mu;
L41c_Lbar_t = L41c_L * mu_ref/mu;

%second segment: from 1 to 2: curved

L12c_index_up = [stringers_root.index_up1 boom_root.index_up(1)];

for i = 1:length(L12c_index_up)-1
    index_seg = [L12c_index_up(i) L12c_index_up(i+1)];
    x = airf_root.XZ_up(index_seg(1) : index_seg(2), 1);
    y = airf_root.XZ_up(index_seg(1) : index_seg(2), 2);
    L=0; 
    A=0;
    for j = 1: length(x)-1
        L =  L + sqrt((x(j) - x(j+1))^2 + (y(j) - y(j+1))^2);
        
        h(j) = abs( (x(j+1) - x(j))*y(j) - x(j)*(y(j+1) - y(j)) )/( (x(j+1) - x(j))^2 + (y(j+1) - y(j))^2 )^(1/2);
        A = A + h(j)*1/2 * sqrt((x(j) - x(j+1))^2 + (y(j) - y(j+1))^2);
    end
    L_up(i) = L;
    L12c_A_up(i) = A;
end

L12c_index_low = [1 stringers_root.index_low1 boom_root.index_low(1)];

for i = 1:length(L12c_index_low)-1
    index_seg = [L12c_index_low(i) L12c_index_low(i+1)];
    x = airf_root.XZ_low(index_seg(1) : index_seg(2), 1);
    y = airf_root.XZ_low(index_seg(1) : index_seg(2), 2);
    L=0; 
    A=0;
    for j = 1: length(x)-1
        L =  L + sqrt((x(j) - x(j+1))^2 + (y(j) - y(j+1))^2);
        
        h(j) = abs( (x(j+1) - x(j))*y(j) - x(j)*(y(j+1) - y(j)) )/( (x(j+1) - x(j))^2 + (y(j+1) - y(j))^2 )^(1/2);
        A = A + h(j)*1/2 * sqrt((x(j) - x(j+1))^2 + (y(j) - y(j+1))^2);
    end
    L_low(i) = L;
    L12c_A_low(i) = A;
end

L12c_l = [flip(L_up) L_low];
L12c_A = [flip(L12c_A_up) L12c_A_low];
L12c_L = sum(L12c_l);
L12c_lbar_t = L12c_l * mu_ref/mu;
L12c_Lbar_t = L12c_L * mu_ref/mu;


%third segment: from 1 to 2: straight line

x = [boom_root.XZ_up(1,1) boom_root.XZ_low(1,1)];
y = [boom_root.XZ_up(2,1) boom_root.XZ_low(2,1)];
L12s_L = sqrt((x(1) - x(2))^2 + (y(1) - y(2))^2);
L12s_Lbar_t = L12s_L * mu_ref/mu;


%fourth segment: from 2 to 3: curve

L23c_index = [boom_root.index_low(1) stringers_root.index_low2 boom_root.index_low(2)];

for i = 1:length(L23c_index)-1
    index_seg = [L23c_index(i) L23c_index(i+1)];
    x = airf_root.XZ_low(index_seg(1) : index_seg(2), 1);
    y = airf_root.XZ_low(index_seg(1) : index_seg(2), 2);
    L=0; 
    A=0;
    for j = 1: length(x)-1
        L =  L + sqrt((x(j) - x(j+1))^2 + (y(j) - y(j+1))^2);
        
        h(j) = abs( (x(j+1) - x(j))*y(j) - x(j)*(y(j+1) - y(j)) )/( (x(j+1) - x(j))^2 + (y(j+1) - y(j))^2 )^(1/2);
        A = A + h(j)*1/2 * sqrt((x(j) - x(j+1))^2 + (y(j) - y(j+1))^2);
    end
    L23c_l(i) = L;
    L23c_A(i) = A;
end

L23c_L = sum(L23c_l);
L23c_lbar_t = L23c_l * mu_ref/mu;
L23c_Lbar_t = L23c_L * mu_ref/mu;


%fifth segment: from 3 to 4: straight

x = [boom_root.XZ_up(1,2) boom_root.XZ_low(1,2)];
y = [boom_root.XZ_up(2,2) boom_root.XZ_low(2,2)];
L34s_L = sqrt((x(1) - x(2))^2 + (y(1) - y(2))^2);
L34s_Lbar_t = L34s_L * mu_ref/mu;


%sixth segment: from 3 to 4: curve

x = airf_root.XZ_up(boom_root.index_up(2) : end, 1);
y = airf_root.XZ_up(boom_root.index_up(2) : end, 2);
L=zeros(length(x),1);
A=0;
for j = 1: length(x)-1
    L(j) =  L(j) + sqrt((x(j) - x(j+1))^2 + (y(j) - y(j+1))^2);
    
        h(j) = abs( (x(j+1) - x(j))*y(j) - x(j)*(y(j+1) - y(j)) )/( (x(j+1) - x(j))^2 + (y(j+1) - y(j))^2 )^(1/2);
        A = A + h(j)*1/2 * sqrt((x(j) - x(j+1))^2 + (y(j) - y(j+1))^2);
end
L_up = L;
L34c_A_up = A;

x = airf_root.XZ_low(boom_root.index_low(2) : end, 1);
y = airf_root.XZ_low(boom_root.index_low(2) : end, 2);
L=zeros(length(x),1);
A=0;
for j = 1: length(x)-1
    L(j) =  L(j) + sqrt((x(j) - x(j+1))^2 + (y(j) - y(j+1))^2);
    
        h(j) = abs( (x(j+1) - x(j))*y(j) - x(j)*(y(j+1) - y(j)) )/( (x(j+1) - x(j))^2 + (y(j+1) - y(j))^2 )^(1/2);
        A = A + h(j)*1/2 * sqrt((x(j) - x(j+1))^2 + (y(j) - y(j+1))^2);
end
L_low = L;
L34c_A_low = A;

L34c_A = [L34c_A_low L34c_A_up];
L34c_l = [L_low; flip(L_up)];
L34c_L = sum(L34c_l);
L34c_Lbar_t = L34c_L * mu_ref/mu;



%Cell 1
Cell_1_L = L12s_L + L12c_L;
Cell_1_Lbar_t = L12s_Lbar_t + L12c_Lbar_t;
Cell_1_area = trapz(airf_root.XZ_up(L12c_index_up, 1),airf_root.XZ_up(L12c_index_up, 2)) - trapz(airf_root.XZ_low(L12c_index_low, 1), airf_root.XZ_low(L12c_index_low, 2));

%Cell 2 
Cell_2_L = L34s_L + L23c_L + L12s_L + L41c_L;
Cell_2_Lbar_t = L34s_Lbar_t + L23c_Lbar_t + L12s_Lbar_t + L41c_Lbar_t;
Cell_2_area = trapz(airf_root.XZ_up(L41c_index(1) : L41c_index(end), 1), airf_root.XZ_up(L41c_index(1) : L41c_index(end), 2)) - trapz(airf_root.XZ_low(L23c_index(1) : L23c_index(end), 1), airf_root.XZ_low(L23c_index(1) : L23c_index(end), 2));

%Cell 3 
Cell_3_L = L34s_L + L34c_L;
Cell_3_Lbar_t = L34s_Lbar_t + L34c_Lbar_t;
Cell_3_area = trapz(airf_root.XZ_up([boom_root.index_up(2) : end], 1), airf_root.XZ_up([boom_root.index_up(2) : end], 2)) - trapz(airf_root.XZ_low([boom_root.index_low(2) : end], 1), airf_root.XZ_low([boom_root.index_low(2) : end], 2));


%% Corrections
syms q_t_1 q_t_2 q_t_3 dtheta_x_t

%Twist rate on Cell 1
eq1 = dtheta_x_t == 1/(2*Cell_1_area * mu_ref) * (q_t_1*Cell_1_Lbar_t - q_t_2*L12s_Lbar_t  + sum(q0_12c .* L12c_lbar_t) - q0_12s*L12s_Lbar_t );

%Twist rate on cell 2
eq2 = dtheta_x_t == 1/(2*Cell_2_area * mu_ref) *(q_t_2*Cell_2_Lbar_t - q_t_1*L12s_Lbar_t - q_t_3*L34s_Lbar_t + q0_12s*L12s_Lbar_t + q0_34s*L34s_Lbar_t + sum(q0_41c.*L41c_lbar_t) + sum(q0_23c.*L23c_lbar_t));

%Twist rate on cell 3
eq3 = dtheta_x_t == 1/(2*Cell_3_area * mu_ref) *(q_t_3*Cell_3_Lbar_t - q_t_2*L34s_Lbar_t - q0_34s*L34s_Lbar_t);


M1 = sum(2*L12c_A.*q0_12c) + sum(2*L23c_A.*q0_23c) + sum(2*L41c_A.*q0_41c) - q0_12s*boom_root.XZ_low(1,1)*L12s_L + q0_34s*boom_root.XZ_low(1,2)*L34s_L;
M2 = 2*Cell_1_area*q_t_1 + 2*Cell_2_area*q_t_2 + 2*Cell_3_area*q_t_3;
M3 = (sum(boom_root.XZ_up(1,:).*boom.P_z_up) + sum(boom_root.XZ_low(1,:).*boom.P_z_up) + sum(stringers_root.XZ_up1(1,:).*stringers.P_z_up1) + sum(stringers_root.XZ_up2(1,:).*stringers.P_z_up2) + sum(stringers_root.XZ_low1(1,:).*stringers.P_z_low1) + sum(stringers_root.XZ_low2(1,:).*stringers.P_z_low2));
M4 = -(sum(boom_root.XZ_up(2,:).*boom.P_x_up) + sum(boom_root.XZ_low(2,:).*boom.P_x_up) + sum(stringers_root.XZ_up1(2,:).*stringers.P_x_up1) + sum(stringers_root.XZ_up2(2,:).*stringers.P_x_up2) + sum(stringers_root.XZ_low1(2,:).*stringers.P_x_low1) + sum(stringers_root.XZ_low2(2,:).*stringers.P_x_low2));

eq4 = M_x == M1 + M2 + M3 + M4;

%Solution
sol = solve([eq1, eq2, eq3, eq4], [q_t_1, q_t_2, q_t_3, dtheta_x_t]);
qt1 = double(sol.q_t_1);
qt2 = double(sol.q_t_2);
qt3 = double(sol.q_t_3);
dthetaxt = double(sol.dtheta_x_t );


%% Shear flow in total:

q_41c = q0_41c + qt2;
q_12c = q0_12c + qt1;
q_12s = q0_12s + qt2 - qt1;
q_23c = q0_23c + qt2;
q_34s = q0_34s + qt2 - qt3;
q_34c = qt3;

q_max = max( [ max(abs(q_41c)) max(abs(q_12c)) max(abs(q_12s)) max(abs(q_23c)) max(abs(q_34s)) max(abs(q_34c))]);

thickness = q_max/tau_max;                   %minimum thickness [m]   



end

