function [] = Maneuver_envelope(in, n, Ue)
%At design altitude: 30 000ft

W = in.MTOW * 9.81; %MTOW in N
rho = density (30000); %[kg/m^3]
rho_0 = 1.225; % [kg/m^3]

[Vc, Vd] = Placard_diagram(); %[m/s] max thrust cruise speed and dive speed
Ve = 0 : 0.1 : Vd; % [m/s] Equivalent airspeed
%stall_line_flaps_up = zeros(1,length(Ve));

%conversion factors to change units
m_ft = 3.28084; %convertion from m to ft
ms_kt = 1.94384; %convertion fromm/s to kt
m2_ft2 = 10.7639104; %from m^2 to ft^2
N_lbf = 0.2248089431; %from N to lb-force
kgm3_slft3 = 0.0019403203; %from kg/m^3 to slug/ft^3
kt_fts = 1.68780986;%from kt to ft/s

%computation of the curves
n_cut = linspace(0, n.max, length(Ve));
stall_line_flaps_up = rho_0 * Ve.^2 * in.S * in.C_L_max /(2*W);
stall_line_flaps_up1 = stall_line_flaps_up;
stall_line_landing = rho_0 * Ve.^2 * in.S * in.C_L_max0 /(2*W);

n_max_line = zeros(1, length(Ve));
n_max_line = n_max_line + 1*n.max;
n_min_line = zeros(1, length(Ve));
n_min_line2 = zeros(1, length(Ve));
n_min_line = n_min_line + 1*n.min;
n_Vd_line = zeros(1,length(Ve));
n_Vd_line = n_Vd_line + 1*Vd;

%cutting the lines to make an envelope
indice1 = find(abs(n_max_line - stall_line_flaps_up)<0.01);
indice2 = find(abs(n_min_line - stall_line_landing)<0.01);

stall_line_flaps_up(indice1(1):end) = [];
Ve1 = Ve(1:indice1(1)-1);
stall_line_landing(indice2:end) = [];
Ve2 = Ve(1:indice2-1);
n_max_line(1:indice1) = [];
Ve3 = Ve(indice1+1:end);
n_min_line(1:indice2) = [];
Ve4 = Ve(indice2+1:end);

%n min line on the right part of the graph
x = [Vd*ms_kt 0];
y = [Vc*ms_kt -1];
m=(y(2)-x(2))/(y(1)-x(1));
p = x(2)-m*x(1);
n_min_line2 =m*Ve*ms_kt +p;

cut_n = find(abs(n_min_line2-n.min)<0.01);
n_min_line2(1:cut_n(1)) = [];
Ve5 = Ve(cut_n(1)+1:end);

cut_nmin = find(abs(Ve5(1)-Ve4)<0.01);
n_min_line(cut_nmin(1)+1:end)= [];
Ve4(cut_nmin(1)+1:end)= [];

%% Plot of the maneuver envelope

%from m/s to kt for x axis
Ve1 = Ve1 * ms_kt;
Ve2 = Ve2 * ms_kt;
Ve3 = Ve3 * ms_kt;
Ve4 = Ve4 * ms_kt;
Ve5 = Ve5 * ms_kt;
Vd = Vd * ms_kt;
Vc = Vc * ms_kt;
Ve = Ve * ms_kt; 
n_Vd_line = n_Vd_line *ms_kt;

Figure1=figure(1); clf; set(Figure1,'defaulttextinterpreter','latex');
hold on
plot(Ve1, stall_line_flaps_up,'linewidth', 2, 'MarkerSize', 11', 'color', '#FF9E00')
plot(Ve2, stall_line_landing,'linewidth', 2, 'MarkerSize', 11', 'color',  '#FF9E00')
plot(Ve3, n_max_line, 'linewidth', 2, 'MarkerSize', 11,'color',  '#FF9E00')
plot(Ve4, n_min_line, 'linewidth', 2, 'MarkerSize', 11,'color',  '#FF9E00')
plot(n_Vd_line, n_cut, 'linewidth', 2, 'MarkerSize', 11,'color', 'black')
plot(Ve5, n_min_line2, 'linewidth', 2, 'MarkerSize', 11,'color', '#FF9E00')

xlabel('Equivalent airspeed [kt]','Fontsize',11)
ylabel('Load factor n [-]','Fontsize',11)
box on
set(gca,'fontsize',11,'fontname','Times', 'LineWidth',0.5);
set(gca,'XMinorTick','off','YMinorTick','off')
set(gca, 'YLim', [-3.5 6]);

%% Gust alleviation factor, using FAR simple rule
%all datas are changed into imperial units to use FAR rule

in.S = in.S*m2_ft2; %from m^2 to ft^2
W = W * N_lbf; %from N to lb-force
in.cw_MAC = in.cw_MAC * m_ft; %from m to ft
rho = rho*kgm3_slft3; %from kg/m^3 to slug/ft^3
rho_0 = rho_0 * kgm3_slft3;
g = 32.174; %gravity in ft/s^2

mu = (2*W)/(rho*in.C_L_alpha_plane*in.cw_MAC*g*in.S); %airplane weight ratio [-]
F = (0.88*mu)/(5.3+mu); %gust alleviation factor 0.8/0.9 [-]

%2 upper lines of the envelope
ng_Vb_line = 1 + ((F * in.C_L_alpha_plane * Ue.Vb * Ve * in.S)/(498*W));
ng_Vb_neg_line = 1 + ((F * in.C_L_alpha_plane * (-1) * Ue.Vb * Ve * in.S)/(498*W));

ng_Vc = 1 + ((F * in.C_L_alpha_plane * Ue.Vc * Vc * in.S)/(498*W));
ng_Vc_neg = 1 + ((F * in.C_L_alpha_plane *(-1) * Ue.Vc * Vc * in.S)/(498*W));

ng_Vd = 1 + ((F * in.C_L_alpha_plane * Ue.Vd * Vd * in.S)/(498*W));
ng_Vd_neg = 1 + ((F * in.C_L_alpha_plane * (-1) * Ue.Vd * Vd * in.S)/(498*W));

ng_Vc_line = 1 + ((F * in.C_L_alpha_plane * Ue.Vc * Ve * in.S)/(498*W));
ng_Vc_neg_line = 1 + ((F * in.C_L_alpha_plane *(-1) * Ue.Vc * Ve * in.S)/(498*W));

ng_Vd_line = 1 + ((F * in.C_L_alpha_plane * Ue.Vd * Ve * in.S)/(498*W));
ng_Vd_neg_line = 1 + ((F * in.C_L_alpha_plane * (-1) * Ue.Vd * Ve * in.S)/(498*W));

indice3 = find(abs(ng_Vb_line - stall_line_flaps_up1)<0.01);
indice3(2:end) = [];
ng_Vb = ng_Vb_line(indice3);
Vb = Ve(indice3); %je pense que c Vb ici et donc on a plus de V5

indice4 = find(abs(ng_Vb_neg_line - stall_line_flaps_up1)<0.01);
indice4(2:end) = [];
ng_Vb_neg = ng_Vb_neg_line(indice4);
Vb_neg = Ve(indice4);

 %{
%cutting the excess lines on the left
indice3 = find(abs(ng_Vb - stall_line_flaps_up1)<0.01);
ng_Vb(1:indice3) = [];
Ve5 = Ve(indice3+1:end);
indice4 = find(abs(ng_Vb_neg - stall_line_flaps_up1)<0.01);
ng_Vb_neg(1:indice4) = [];
Ve6 = Ve(indice4+1:end);

%cutting the excess lines on the right
indice5 = find(abs(Ve5-Vc)<0.1);
ng_Vb(indice5:end) = [];
Ve5 = Ve5(1:indice5-1);
indice6 = find(abs(Ve6-Vc)<0.1);
ng_Vb_neg(indice6:end) = [];
Ve6 = Ve6(1:indice6-1);
%}

stall_line_flaps_up_gust = stall_line_flaps_up(indice4:indice3);
Ve_stall_gust = Ve(indice4:indice3);

Ve_cruise_index = find(abs(stall_line_flaps_up1 - 1)<0.01);
Ve_cruise = Ve(Ve_cruise_index(1));
Vs1 = Ve_cruise;

%Vb = Ve_stall_gust(end);
%% plot the additionnal curves for the gust envelope
hold on
plot([Vb Vc], [ng_Vb ng_Vc],'linewidth', 1.5, 'MarkerSize', 11', 'color', '#00707F')
plot([Vc Vd], [ng_Vc ng_Vd],'linewidth', 1.5, 'MarkerSize', 11', 'color',  '#00707F')
plot([Vb_neg Vc], [ng_Vb_neg ng_Vc_neg],'linewidth', 1.5, 'MarkerSize', 11', 'color', '#00707F')
plot([Vc Vd], [ng_Vc_neg ng_Vd_neg],'linewidth', 1.5, 'MarkerSize', 11', 'color',  '#00707F')
plot(Ve_stall_gust, stall_line_flaps_up_gust,'linewidth', 1.5, 'MarkerSize', 11', 'color', '#00707F','LineStyle',':')

plot(Ve, ng_Vb_line,'linewidth', 1, 'MarkerSize', 11', 'color',  'black','LineStyle',':')
plot(Ve, ng_Vc_line,'linewidth', 1, 'MarkerSize', 11', 'color',  'black','LineStyle',':')
plot(Ve, ng_Vd_line,'linewidth', 1, 'MarkerSize', 11', 'color',  'black','LineStyle',':')
plot(Ve, ng_Vc_neg_line,'linewidth', 1, 'MarkerSize', 11', 'color',  'black','LineStyle',':')
plot(Ve, ng_Vd_neg_line,'linewidth', 1, 'MarkerSize', 11', 'color',  'black','LineStyle',':')


%plot of the velocity legend
plot([Vd Vd], [n.max -2.5],'linewidth', 1, 'MarkerSize', 11', 'color',  'black','LineStyle',':')
plot([Vc Vc], [ng_Vc -2.5],'linewidth', 1, 'MarkerSize', 11', 'color',  'black','LineStyle',':')
plot([Vb Vb], [stall_line_flaps_up_gust(end) -2.5],'linewidth', 1, 'MarkerSize', 11', 'color',  'black','LineStyle',':')
plot([Vs1 Vd], [1 1],'linewidth', 1, 'MarkerSize', 11', 'color',  'black','LineStyle','--')
plot([Vs1 Vs1], [1 -2.5],'linewidth', 1, 'MarkerSize', 11', 'color',  'black','LineStyle',':')
plot([Ve1(end) Ve1(end)], [stall_line_flaps_up(end) -2.5],'linewidth', 1, 'MarkerSize', 11', 'color',  'black','LineStyle',':')

plot([0 Ve(end)], [n.max*1.5 n.max*1.5],'linewidth', 1, 'MarkerSize', 11', 'color',  'black','LineStyle','--')
plot([0 Ve(end)], [n.max n.max],'linewidth', 1, 'MarkerSize', 11', 'color',  'black','LineStyle','--')
plot([0 Ve(end)], [n.min n.min],'linewidth', 1, 'MarkerSize', 11', 'color',  'black','LineStyle','--')

text(250,1.23,'Cruise','interpreter','latex','Fontsize',11)
text(Vd-20,-2,'Vd','interpreter','latex','Fontsize',11)
text(Vc-20,-2,'Vc','interpreter','latex','Fontsize',11)
text(Vs1-20,-2,'Vs1','interpreter','latex','Fontsize',11)
text(Vb-20,-2,'Vb','interpreter','latex','Fontsize',11)
text(Ve1(end)-20,-2,'Va','interpreter','latex','Fontsize',11)

text(20,(n.max*1.5)+0.2,'$n_{ultimate}$','interpreter','latex','Fontsize',11)
text(20,(n.max)+0.2,'$n_{max}$','interpreter','latex','Fontsize',11)
text(20,n.min+0.2,'$n_{min}$','interpreter','latex','Fontsize',11)

hgexport(Figure1,'Maneuver_envelope');

%% FAR requirements slide 20
if Vb < (Vs1*ng_Vc^(1/2))
    disp('ok 1')
else
    disp('requirement 1 is not fulfilled')
end

Ue = [Ue.Vc Ue.Vb Ue.Vd];

if Vc > Vb + 1.32*Ue
    disp('ok 2')
else
    disp('requirement 2 is not fulfilled')
end

%Vc = Vc* kt_fts;%from kt to ft/s
test3 = Vs1 * sqrt(1+(rho_0 * Vc * in.S * F * in.C_L_alpha_plane * Ue)/(2*W));
if Vb > test3
    disp('ok 3')
else
    disp('requirement 3 is not fulfilled')
end

end

