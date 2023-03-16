function [V_c, V_d] = Placard_diagram()
    %PLACARD_DIAGRAM returns the equivalent cruise and dive speed.
    %It also plots the placard diagram of the airplane.
    % -----
    %
    % Syntax:
    %   [V_c, V_d] = Placard_diagram(alt) returns an array with
    %   the equivalent cruise speed [m/s] and the equivalent dive speed [m/s].
    %
    % Outputs:
    %   V_c: the equivalent cruise speed of the aircraft[m/s]
    %   V_d: the equivalent dive speed of the aircraft [m/s] 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data
T_max = 3.36*10^3;                  %Maximum thrust provided by the engine[N]
S_w = 119.5/10.764;                 %Surface of the wing[m²]
AR_w = 7;                           %Aspect ratio of the wing
C_d0 = 0.017;                      %Drag independant of the lift
e = 0.8;                            %Oswald efficiency factor
W = 4228*9.18;                      %Weight [N]

rho_0 = 1.225;                      %air density at sea level [kg/m^3]
rho_design = 0.4594;                %air density at design altitude [kg/m^3] %cfr table florian

speed_range = 0:0.1:300;                                                        %[m/s]
design_alt = 30000*0.0003048*ones(size(speed_range));                           %[km]
turbulence_zone = [10.6*ones(size(speed_range)) ;7.6*ones(size(speed_range))];  %[km]
stratosphere_limit = 11*ones(size(speed_range));                                %[km]

%%Cruise
%%true_speed at design altitude
speed_design = 240; %[m/s]
error = 1;
tol = 10^-3;

while(error > tol)
    speed_design_old = speed_design;
    C_L = W/(0.5 * rho_design * speed_design^2 * S_w);
    C_D = C_d0 + (C_L^2)/(e*pi*AR_w);
    speed_design = sqrt(T_max/(0.5*rho_design*S_w*C_D));
    error = abs((speed_design_old - speed_design)/speed_design_old);
end

a_design = speed(30000,1);                                 %speed of sound at design altitude[m/s]
M_c = speed_design/a_design;                               %mach at design altitude

speed_strat = speed(stratosphere_limit(1)/0.0003048, M_c); %equivalent velocity at stratosphere altitude [m/s]
speed_sea = speed_design*sqrt(rho_design/rho_0);           %equivalent velocity at sea level [m/s]


speed_C = [speed_sea speed_design speed_strat speed_strat];
altitude_C = [0 design_alt(1) stratosphere_limit(1) stratosphere_limit(1)+3];
altitude_C = altitude_C/0.0003048;                          %[ft]
speed_C = speed_C*1.9438444924;                             %[kts]

%%Dive
M_d = 1.07*M_c;                                              %Mach at dive[-]
speed_design_D = speed(design_alt(1)/0.0003048, M_d);        %dive speed at design altitude [m/s]
speed_strat_D = speed(stratosphere_limit(1)/0.0003048, M_d); %dive speed at stratosphere altitude [m/s]

%equivalent speed = min(1.15*V_C at sea level and V(M_d))
find_speed = speed(0, M_d);
find_speed_2 = 1.15*speed_sea;
find_alt = 0;
while(find_speed-find_speed_2) > 0.001
find_alt = find_alt + 0.1;
find_speed = speed(find_alt/0.0003048, M_d);
find_rho = density(find_alt/0.0003048);
find_speed_2 = 1.15*speed_sea/sqrt(find_rho/rho_0); 
end

speed_D = [1.15*speed_sea find_speed_2 speed_design_D speed_strat_D speed_strat_D];
altitude_D = [0 find_alt design_alt(1) stratosphere_limit(1) stratosphere_limit(1)+3];
altitude_D = altitude_D/0.0003048;  %[ft]
speed_D = speed_D*1.9438444924;     %[kts]

%% Return the equivalent velocity [m/s]
 V_c = speed_sea; 
 V_d = 1.15*speed_sea;

%% Plot

Figure1=figure(1); clf; set(Figure1,'defaulttextinterpreter','latex'); 
hold on

plot(speed_C,altitude_C,'linewidth', 2, 'MarkerSize', 15', 'color', '#00707F')
plot(speed_D,altitude_D,'linewidth', 2, 'MarkerSize', 15', 'color',  '#FF9E00')
plot(speed_range*1.9438444924,design_alt/0.0003048,'linewidth', 2, 'MarkerSize', 15,'color', 'k','LineStyle','--') 
plot(speed_range*1.9438444924,stratosphere_limit/0.0003048,'linewidth', 2, 'MarkerSize', 15,'color', 'k','LineStyle','-.')
plot(speed_range*1.9438444924,turbulence_zone/0.0003048,'linewidth', 2, 'MarkerSize', 15,'color', 'k','LineStyle',':')

text(0,30000+800,'Design altitude','interpreter','latex','Fontsize',24)
text(0,turbulence_zone(2, 1)/0.0003048+800,'Turbulence zone','interpreter','latex','Fontsize',24)
text(0,stratosphere_limit(1)/0.0003048+800,'Stratosphere limit','interpreter','latex','Fontsize',24)

xlabel('True airspeed [kts]','Fontsize',24)
ylabel('Altitude [ft]','Fontsize',24)
box on
legend('Cruise speed', 'Dive speed','interpreter','latex','Fontsize',24, 'Location', 'southwest');
set(gca,'fontsize',24,'fontname','Times', 'LineWidth',0.5);
set(gca,'XMinorTick','off','YMinorTick','off')
hgexport(Figure1,'placard_diagram');
end