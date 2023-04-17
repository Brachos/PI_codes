feet = 3.28084; % m to ft
h = (0:10:30000)/feet; % [m] Altitude

A = zeros(1,length(h));
Z = zeros(1,length(h));
X = zeros(1,length(h));
P = zeros(1,length(h)); % [Pa] Pressure(altitude)

T0 = 13580; % [N] Static thrust
BPR = 4.1;
G = 1.1;
Pb = 101325; % [Pa] Pressure at see level
Tb = 273.15 + 15; % [K] Temperature at see level
hb = 11000; % [m] Height at the bottom of atmospheric layer
Lb = -0.0065; % [K/m]
R = 8.3142; % [N.m/mol.K]
g = 9.80665; % [m/s^2]
MM = 0.0289644; % [kg/mol] Molar mass of air
[~,loc] = min(abs(h-hb));

for i = 1 : length(h)
    if i <= loc
        P(i) = Pb*(1 + Lb/Tb*(h(i)))^(-g*MM/(R*Lb));
    else
        P(i) = P(loc)*exp(-g*MM*(h(i)-h(loc))/(R*(Tb-71.5)));
    end
    A(i) = -0.4327*(P(i)/Pb)^2 + 1.3855*(P(i)/Pb) + 0.0427;
    X(i) = 0.1377*(P(i)/Pb)^2 - 0.4374*(P(i)/Pb) + 1.3003;
    Z(i) = 0.9106*(P(i)/Pb)^2 - 1.7736*(P(i)/Pb) + 1.8697;

end

M = [0 0.1 0.2 0.5 0.7];
T = zeros(length(M),length(h));

figure
hold on
for i = 1 : length(M)
    for j = 1 : length(h)
        T(i,j) = T0*(A(j) - 0.377*(1+BPR)*Z(j)/sqrt(G*(1+0.82*BPR))*P(j)/Pb*M(i) + (0.23+0.19*sqrt(BPR))*X(j)*P(j)/Pb*M(i)^2);
    end
    plot(T(i,:)/1000,h*feet)
end
hold off
legend('M = 0','M = 0.1','M = 0.2','M = 0.5','M = 0.7')
ylabel('Altitude [ft]')
xlabel('Thrust [kN]')
