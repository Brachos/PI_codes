function [rho] = density (altitude)
    %DENSITY returns the air density at altitude.
    % -----
    %
    % Syntax:
    %   [speed] = speed(altitude,mach) returns the air density [kg/m³] at altitude [ft].
    %
    % Inputs:
    %   altitude: the altitude at which the speed and the air density must
    %   be computed [ft]
    %
    % Outputs:
    %   rho: the air density at altitude [kg/m³]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
T0 = 15 + 273.15;                       %temperature at sea level
delta_T = altitude / 1000 * 2;          %decrease of 2 degrees every 1000ft
T = T0 - delta_T;                       %temperature at altitude [K]
R  = 8314.32;
P0 = 101325;                            %air pressure at sea level [Pa]
g  = 9.81;                              %gravity [m/s²]
z  = altitude*0.3048;                   %[m]
Cp = 1006;
P  = P0*exp(-7*g*z/(2*Cp*T0));          %air pressure at altitude [Pa]
rho = P*28.9644/(R*T);                  %air density at altitude (perfect gases)
end