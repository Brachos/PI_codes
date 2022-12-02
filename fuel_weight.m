function [fuel_mass] = fuel_weight()
%%Quantity of fuel in kg

ingress_range = 850 * 1.852; % [nm] to [km] at 30 000 ft % KPP02 O
ingress_speed = speed(30000,0.7); %[mach] to [km/h] at 30 000ft %KPP03 T
loiter_time = 5;% [hours] at 30 000 ft and 0.7 mach %KPP4 O
[loiter1_speed,rho1] = speed(30000,0.7);
egress_range = 2000 * 1.852; % [nm] to [km] at 30 000 ft % KPP05 O
[egress_speed,rho2] = speed(30000, 0.7); %[mach] to [km/h] at 30 000ft %KPP06 T
loiter_time_landing = 0.75; %[hours] %KPP10 O
[loiter2_speed,rho3] = speed(30000,0.7);
TSFC = 0.486 *0.45 / 4.49; %[lb/lbf/h] to [kg/h/N]
net_thrust = 3600;%[N]

%estimated flying time 
ingress_time = ingress_range / ingress_speed; %[hours] 
egress_time = egress_range / egress_speed; %[hours]
total_flying_time = ingress_time + egress_time + loiter_time + loiter_time_landing;

%reserves

%TSFC = Thrust specific Fuel Consumption
fuel_mass_flow_rate = net_thrust * TSFC; % [N * kg/hr/N] = [kg/hr]
fuel_mass = fuel_mass_flow_rate * total_flying_time; %[kg]
%on a 2000 kg de fuel, et 4500 rempli

%range
loiter1_range = loiter_time * loiter1_speed; %
loiter2_range = loiter_time_landing * loiter2_speed;

range = ingress_range + egress_range + loiter1_range + loiter2_range;
range2 = total_flying_time * ingress_speed; 
end
