function [HE,HT,Hmfg,N_eng,t_ac,CPI,Ceng,Cdev,CFT,Ctool,CMFG,Cqc,Cmat,Ccert,Cpp,cost_per_aircraft,Cstor,Cins,Cinsp,Cfuel,Cap] = Cost_analysis(We,Drag,show)
%% Data required
%WE Aircraft empty weight [lb]
We = We*2,205;
V_max = 492.9; %Maximum air speed in KTAS
%V_max = 250; %Maximum air speed in m/S
T = Drag*0.2248; % Cruise thrust in lbf
%Demander mathias
speed_cruise = 227,44176; %feet/s
%% Assumtion 
N = 100:1:1000; % Number of aircraft produced in 5 years -> 8 aircraft produced per day
%%
pt = 14;
F_EXP = 0.8:0.05:0.95;
Number_unit_produced =1:1:600;
QDF = zeros(length(F_EXP),length(Number_unit_produced)); %Quantity discount factor vector
for i=1:length(F_EXP)
    for j=1:length(Number_unit_produced)
        QDF(i,j) = (F_EXP(i))^(1.4427*log(Number_unit_produced(j)));
    end
end
%% CPI computation
year = 2013:1:2028;
CPI = [1.015 1.016 1.001 1.013 1.021 1.024 1.018 1.012 1.047 1.086 1.0237 1.0237 1.0237 1.0237 1.0237 1.0237];
CPI_evolution = [];
CPI_evolution(1) = 1; %Correspond à 2012
for i=2:length(CPI)
    CPI_evolution(i) = CPI_evolution(i-1)*CPI(i);
end
if show == 1 
    figure; 
    b = plot(year,CPI_evolution,'LineWidth',2);
    xlabel('Time [years]','Interpreter','latex','FontSize',pt);
    ylabel('CPI [-]','Interpreter','latex','FontSize',pt);
    grid on
    box on
    pbaspect([1.5 1 1])
    saveas(b,'.\CPI','epsc')
end
%% Plot of the quantity discount factor
if show == 1
    d=1;
    figure; 
    a = plot(Number_unit_produced,QDF(1,:),'LineWidth',2);
    hold on 
    plot(Number_unit_produced,QDF(2,:),'LineWidth',2)
    plot(Number_unit_produced,QDF(3,:),'LineWidth',2)
    plot(Number_unit_produced,QDF(4,:),'LineWidth',2)
    xlabel('Number of unit produced, N [-]','Interpreter','latex','FontSize',pt);
    ylabel('Quantity discount factor [-]','Interpreter','latex','FontSize',pt);
    grid on
    box on
    h = legend('80$\%$','85$\%$','90$\%$','95$\%$','Reference','interpreter','latex','Fontsize',pt)
    pbaspect([1.5 1 1])
    saveas(a,'.\Quantity_discount_factor','epsc')
end

%% General aviation computation 
for i=1:length(N)
    [HE,HT,Hmfg,N_eng,t_ac,CPI,Ceng,Cdev,CFT,Ctool,CMFG,Cqc,Cmat,Ccert,Cpp,cost_per_aircraft] = general_aviation(We,V_max,400,T);
    cost_per_unit(i) = cost_per_aircraft;
end
if show == 1
    figure;
    b = plot(N,cost_per_unit,'LineWidth',2);
    xlabel('Number of units produced N [-]','Interpreter','latex','FontSize',pt)
    ylabel('Unit selling price [$]','Interpreter','latex','FontSize',pt)
end 
%% Break even analysis 

%% Operational costs
F1 = -0.15; %Maintenance performed by owner 
F2 = 0.02; %Difficult engine access
F3 = 0.02; %Retreactable landing gear
F4 = 0.02; %VFR radios are installed
F5 = 0.04; %IFR radio installed 
F6 = 0.01; %Integral fuel tank installed
F7 = 0; %Simple flap system
F8 = 0; % 14 CFR Part 23

FMF = 0.3 +F1+F2+F3+F4+F5+F6+F7+F8;
Rap = 60; %dollar
Rstor = 250; %Storage price [$]
Cstor = 12 * Rstor*CPI; %Storage cost (annual) [$]

Rfuel = 2.68; % Price of kerozene per gallon (voir source) [$]
SFC = 0.545; %Specific fuel consumption at cruise
BHP_cruise = T*speed_cruise;
Qflght = 300; %Flight hours per year
Cap = FMF*Rap*Qflght;
Cfuel = BHP_cruise*SFC*Qflght*Rfuel;

Cac = cost_per_aircraft; %Total value assured [$]
Cins = 500 + 0.015*Cac; % Annual insurance cost [$/year]
Cinsp = 500*CPI; %Anual inspection cost
Npp = 1; %Number of engine
Cover = 5*Npp*Qflght; %Engine overhaul fund





end 