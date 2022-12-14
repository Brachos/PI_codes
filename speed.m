function [speed,rho] = speed(altitude,mach)

delta_T = altitude / 1000 * 2; %diminution de 2 degrés par 1000 ft
K = 273.15;
T = 15 - delta_T + K;
R_s = 287.058; %[J/kg/K] constante spécifique de l'air sec
gamma = 1.4; %hypothèse des gaz parfaits, coeff de Laplace
speed = mach * sqrt(gamma * T * R_s);%vitesse [m/s]

R  = 8314.32;
P0 = 101325;
T0 = 15+K;
g  = 9.81;
z  = altitude*0.3048;
Cp = 1006;
P  = P0*exp(-7*g*z/(2*Cp*T0));

rho = P*28.9644/(R*T); % Gaz parfait

end

