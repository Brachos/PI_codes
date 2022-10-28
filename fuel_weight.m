%%Quantity of fuel in kg

ingress_range = 650 * 1.852; % [nm] to [km] at 30 000 ft % KPP02
ingress_speed = speed(30000,0.7); %[mach] to [km/h] at 30 000ft %KPP03
loiter_time = 3.5;% [hours] at 30 000 ft and 0.7 mach %KPP4
egress_range = 1650 * 1.852; % [nm] to [km] at 30 000 ft % KPP05
egress_speed = speed(30000, 0.5); %[mach] to [km/h] at 30 000ft %KPP06
loiter_time_landing = 0.5; %[hours]
TSFC = 0.7 *0.45 / 4.49; %[kg/h/N]
net_thrust = 3600;%[N]

%estimated flying time 
ingress_time = ingress_range / ingress_speed; %[hours] 
egress_time = egress_range / egress_speed; %[hours]
total_flying_time = ingress_time + egress_time + loiter_time + loiter_time_landing;

%reserves

%TSFC = Thrust specific Fuel Consumption
fuel_mass_flow_rate = net_thrust * TSFC; % [N * kg/hr/N] = [kg/hr]
fuel_mass = fuel_mass_flow_rate * total_flying_time %[kg]
%on a 2000 kg de fuel, et 4500 rempli

