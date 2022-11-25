clear
M_cruise = 0.7; %[-]
M_c = 1.06*M_cruise; %[-]
M_d = 1.07*M_c; %[-]
rho_0 = 1.225;
altitude = 0:0.1:15; %[km]
speed_range = 0:0.1:300; %[m/s]
design_alt = 30000*0.0003048*ones(size(speed_range)); %[km]
turbulence_zone = [10.6*ones(size(speed_range)) ;7.6*ones(size(speed_range))]; %[km]
stratosphere_limit = 11*ones(size(speed_range)); %[km]

%%Cruise
[speed_true_design, rho_design] = speed(design_alt(1)/0.0003048, M_c);
[speed_true_strat, rho_strat] = speed(stratosphere_limit(1)/0.0003048, M_c);
speed_eq = speed_true_design*sqrt(rho_design/rho_0); %at 0[m] v_true = v_eq;

speed_C = [speed_eq speed_true_design speed_true_strat speed_true_strat];
altitude_C = [0 design_alt(1) stratosphere_limit(1) stratosphere_limit(1)+3];

%%Dive
[speed_true_design_D, rho_design_D] = speed(design_alt(1)/0.0003048, M_d);
[speed_true_strat_D, rho_strat_D] = speed(stratosphere_limit(1)/0.0003048, M_d);

find_speed = speed(0, M_d);
find_speed_2 = 1.15*speed_eq;
find_alt = 0;
while(find_speed-find_speed_2) > 0.001
find_alt = find_alt + 0.1;
[find_speed, find_rho] = speed(find_alt/0.0003048, M_d);
find_speed_2 = 1.15*speed_eq/sqrt(find_rho/rho_0); 
end

speed_D = [1.15*speed_eq find_speed_2 speed_true_design_D speed_true_strat_D speed_true_strat_D];
altitude_D = [0 find_alt design_alt(1) stratosphere_limit(1) stratosphere_limit(1)+3];

%% Plot
Figure1=figure(1); clf; set(Figure1,'defaulttextinterpreter','latex'); 
hold on
plot(speed_range,design_alt,'linewidth', 2, 'MarkerSize', 15,'color', 'r') 
plot(speed_range,stratosphere_limit,'linewidth', 2, 'MarkerSize', 15,'color', 'g')
plot(speed_range,turbulence_zone,'linewidth', 2, 'MarkerSize', 15,'color', 'b')
plot(speed_C,altitude_C,'linewidth', 2, 'MarkerSize', 15')
plot(speed_D,altitude_D,'linewidth', 2, 'MarkerSize', 15')

xlabel('True airspeed [m/s]')
ylabel('Altitude [km]')
box on
legend('Design altitude','Stratosphere limit','Turbulence Zone','interpreter','latex','Fontsize',30, 'Location', 'northwest');
set(gca,'fontsize',24,'fontname','Times', 'LineWidth',0.5);
set(gca,'XMinorTick','off','YMinorTick','off')
%hgexport(Figure1,'placard_diagram');