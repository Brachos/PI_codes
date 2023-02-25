function [speed,rho] = VelocityCruise(altitude,mach)

%% Parameter defintions
gamma = 1.4;                          % Adiabatic air constant
R  = 287.05;                          % Gas constant in [J/kg.K].
P0 = 101325;                          % Pressure at sea level in [N/m^2].
T0 = 288.15;                          % Temperature at sea level in [K].
g  = 9.80065;                         % Gravity value at sea level in [m/s^2].                               
h  = altitude*0.3048;                 % Ceiling in [m].
k1 = -0.0065;                         % Temperature lapse in [K/m].
%Cp = 1006; 

%% Computations

% h = input('Enter the cruise altitude: ');

T = T0 + k1*h;                                      % Temperature value at cruise in [K].
P  = P0*exp(T/T0)^-(g/(R*k1));                      % pressure value at cruise in [N/m^2]. 

% Velocity and density evaluation at cruise

speed = mach * sqrt(gamma * T * R);                 % velocity at cruise in [m/s].
rho = P/(R*T);                                      % density at cruise in in [kg/m^3].

end