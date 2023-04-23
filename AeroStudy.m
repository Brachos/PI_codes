%% Code for the aerodynamic study 
% Author = Charles Jacquet
%% grF computation
percentage = 0.010345; % wing surface
grF = [1.1 1.2 1.3]; %Growth ratio 

chord_0 = 1.64;
chord_1 = 0.492;
MAC = 1.1714;
meshPerChord = 1/percentage;
boxsize = 70*chord_0;

a = chord_0/meshPerChord;
beta = boxsize*(grF-1)/(2*a);
n = log10(beta)./log10(grF) - 1;

msF = a*grF.^n

%% Convergence results load
Convergence_domain = readtable("Convergence_domain.csv");
Convergence_growth_1_1 = readtable("Convergence_gr1_1.csv");
Convergence_growth_1_2 = readtable("Convergence_gr1_2.csv");
Convergence_growth_1_3 = readtable("Convergence_gr1_3.csv");
pt = 14;
%% Plot convergence domain length = height

figure; 
a = semilogy(Convergence_domain.BoxSizeLength(:,1),Convergence_domain.Cl(:,1),'LineWidth',2);
xlabel('Domain size length = Domain size height [m]','Interpreter','latex','FontSize',pt);
ylabel('Cl [-]','Interpreter','latex','FontSize',pt);
grid on
box on
pbaspect([1.5 1 1])
saveas(a,'.\Convergence_domain','epsc')
%% Plot convergence domain width
figure; 
d = semilogy(Convergence_domain.BoxSizeWidth(:,1),Convergence_domain.Cl(:,1),'LineWidth',2);
xlabel('Domain width [m]','Interpreter','latex','FontSize',pt);
ylabel('Cl [-]','Interpreter','latex','FontSize',pt);
grid on
box on
pbaspect([1.5 1 1])
saveas(d,'.\Convergence_domain_width','epsc')

%% Plot error domain size
figure; 
c = plot(Convergence_domain.BoxSizeLength(:,1),Convergence_domain.Relative_error_percent(:,1),'LineWidth',2);
xlabel('Domain size [m]','Interpreter','latex','FontSize',pt);
ylabel('Relative Error [%]','Interpreter','latex','FontSize',pt);
grid on
box on
pbaspect([1.5 1 1])
saveas(c,'.\Error_domain','epsc')

%% Plot convergence growth ratio 
figure; 
b = semilogy(Convergence_growth_1_1.Number_element(:,1),Convergence_growth_1_1.Cl(:,1),'LineWidth',2);
hold on 
semilogy(Convergence_growth_1_2.Number_element(:,1),Convergence_growth_1_2.Cl(:,1),'LineWidth',2);
semilogy(Convergence_growth_1_3.Number_element(:,1),Convergence_growth_1_3.Cl(:,1),'LineWidth',2);
xlabel('Number of element [-]','Interpreter','latex','FontSize',pt);
ylabel('Cl [-]','Interpreter','latex','FontSize',pt);
grid on
box on
h = legend('gr = 1.1','gr = 1.2','gr = 1.3','interpreter','latex','location','southeast','Fontsize',pt)
pbaspect([1.5 1 1])
saveas(c,'.\growth_ratio','epsc')

%% Convergence taille des elements
Number_element_1_8 = [Convergence_growth_1_3.Number_element(7,1),Convergence_growth_1_2.Number_element(7,1),Convergence_growth_1_1.Number_element(7,1)];
Number_element_1_2 = [Convergence_growth_1_3.Number_element(8,1),Convergence_growth_1_2.Number_element(8,1),Convergence_growth_1_1.Number_element(8,1)];
Number_element_1_03 = [Convergence_growth_1_3.Number_element(9,1),Convergence_growth_1_2.Number_element(9,1),Convergence_growth_1_1.Number_element(9,1)];
Cl_1_8 = [Convergence_growth_1_3.Cl(7,1),Convergence_growth_1_2.Cl(7,1),Convergence_growth_1_1.Cl(7,1)];
Cl_1_2 = [Convergence_growth_1_3.Cl(8,1),Convergence_growth_1_2.Cl(8,1),Convergence_growth_1_1.Cl(8,1)];
Cl_1_03 = [Convergence_growth_1_3.Cl(9,1),Convergence_growth_1_2.Cl(9,1),Convergence_growth_1_1.Cl(9,1)];
figure; 
e=plot(Number_element_1_8,Cl_1_8,'LineWidth',2)
hold on 
plot(Number_element_1_2,Cl_1_2,'LineWidth',2)
plot(Number_element_1_03,Cl_1_03,'LineWidth',2)
xlabel('Number of element [-]','Interpreter','latex','FontSize',pt);
ylabel('Cl [-]','Interpreter','latex','FontSize',pt);
grid on
box on
h = legend('1.8 $\%$','1.2 $\%$','1.03 $\%$','interpreter','latex','location','southeast','Fontsize',pt)
pbaspect([1.5 1 1])
saveas(e,'.\conv_element_size','epsc')

%% Cp_curve
chord = linspace(0.1,3.71,25);
cp_curve_1_1 = load('Converged_mesh_slice_10_MAC_growth_1_1.dat');
% x, y, z, x/c, Cp
figure; 
w = plot(cp_curve_1_1(:,4),cp_curve_1_1(:,5))
set(gca, 'YDir','reverse')
%% Compute of aero coefs
Aero_coefs = readtable("AEROSTUDYCOEFS.csv");
prop_lift = 0.8004;
MTOW = 3605.26;
Mach = 0.7;
Altitude = 30000; % [ft]
Mass = 1000;
AR = 7;
AOA = -1 : 0.1 : 1;
CL_vector = zeros(length(AOA),1);
CD_vector = zeros(length(AOA),1);
deriv = zeros(length(AOA),1);

for i = 1 : length(AOA)
    [~,~,~,~,CL_vector(i),CD_vector(i),~,~,~,~,~,~,~,~,~,~,~,~,~,~] = wing(Mach,Altitude,Mass,AOA(i)+2);
    deriv(i) = 0.5/CL_vector(i)*pi*0.8*AR;
    if abs(deriv(i)*CD_vector(i)-CL_vector(i)) < 0.005
        CL_opt = CL_vector(i);
        num = i;
    end
end
%%
%plot AOA vs Cl 
figure; 
k = plot(Aero_coefs.AOA(:,1),Aero_coefs.CL(:,1),'Linewidth',2)
hold on 
plot(AOA,CL_vector,'LineWidth',2)
xlabel('AoA [$^{\circ}$]','Interpreter','latex','FontSize',pt);
ylabel('$C_L$ [-]','Interpreter','latex','FontSize',pt);
grid on
box on
h = legend('DARTFlo','Emperical','interpreter','latex','location','southeast','Fontsize',pt)
pbaspect([1.5 1 1])
saveas(k,'.\AOAvsCL','epsc')
% Plot AOA vs CD
figure; 
l = plot(Aero_coefs.AOA(:,1),Aero_coefs.CD(:,1),'LineWidth',2)
xlabel('AoA [$^{\circ}$]','Interpreter','latex','FontSize',pt);
ylabel('$C_D$ [-]','Interpreter','latex','FontSize',pt);
grid on
box on
pbaspect([1.5 1 1])
saveas(l,'.\AOAvsCD','epsc')
% Plot AOA vs Cm 
figure; 
m = plot(Aero_coefs.AOA(:,1),Aero_coefs.cm(:,1),'LineWidth',2)
xlabel('AoA [$^{\circ}$]','Interpreter','latex','FontSize',pt);
ylabel('$C_{MAC}$ [-]','Interpreter','latex','FontSize',pt);
grid on
box on
pbaspect([1.5 1 1])
saveas(m,'.\AOAvsCm','epsc')