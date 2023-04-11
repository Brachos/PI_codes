function [speed] = speed(altitude,mach)
    %SPEED returns the speed corresponding to mach at altitude.
    % -----
    %
    % Syntax:
    %   [speed] = speed(altitude,mach) returns an array with
    %   the speed [m/s] at altitude [ft].
    %
    % Inputs:
    %   altitude: the altitude at which the speed and the air density must
    %   be computed [ft]
    %   mach: the mach of the aircraft
    %
    % Outputs:
    %   speed: the speed of the aircraft flying at mach and at altitude
    %   [m/s]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
T0 = 15 + 273.15;                       %temperature at sea level

%speed computation
delta_T = altitude / 1000 * 2;          %decrease of 2 degrees every 1000ft
T = T0 - delta_T;                       %temperature at altitude [K]
R_s = 287.058;                          %dry air specific constant[J/kg/K] 
gamma = 1.4;                            %Laplace coefficient (perfect gases)
speed = mach * sqrt(gamma * T * R_s);   %speed at the design altitude[m/s]


end