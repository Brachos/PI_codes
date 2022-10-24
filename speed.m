function [speed] = speed(altitude,mach)

delta_T = altitude / 1000 * 2; %diminution de 2 degrés par 1000 ft
T = 15 - delta_T + 273.15;
R_s = 287.058; %[J/kg/K] constante spécifique de l'air sec
gamma = 5/3; %hypothèse des gaz parfaits, coeff de Laplace

speed = mach * sqrt(gamma * T * R_s);%vitesse [m/s]
speed = speed * 3.6; %[km/h]

end

