%% Script that plot CL vs AOA and Drag Polar

Mach = 0.7;
Altitude = 30000; % [ft]
Mass = 1000;
AR = 7;
AOA = -10 : 0.1 : 10;
CL_vector = zeros(length(AOA),1);
CD_vector = zeros(length(AOA),1);
deriv = zeros(length(AOA),1);

for i = 1 : length(AOA)
    [~,~,~,~,CL_vector(i),CD_vector(i),~,~,~,~,~,~,~,~,~,~,~,~,~,~] = wing(Mach,Altitude,Mass,AOA(i));
    deriv(i) = 0.5/CL_vector(i)*pi*0.8*AR;
    if abs(deriv(i)*CD_vector(i)-CL_vector(i)) < 0.005
        CL_opt = CL_vector(i);
        num = i;
    end
end

figure
plot(CD_vector,CL_vector)
hold on
plot([0;CD_vector],deriv(135)*[0;CD_vector])
plot(CD_vector(135),CL_vector(135),'xr')
xlabel('CD [-]')
ylabel('CL [-]')

figure 
plot(AOA,CL_vector)
xlabel('AOA [deg]')
ylabel('CL [-]')
