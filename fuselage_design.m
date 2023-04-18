function [Deq_val,a_val,b_val,L_f_val,V_f]=fuselage_design(MTOW,Vw_fuel, net_thrust)
%thickness of the fuselage
t_fuselage=0.05; %[m]

%volume of the payload
V_sensors=0.5;%[m^3]
V_subsystems=0.4;%[m^3]
V_mission=0.15;%[m^3]
L_engine=1.58;
V_engine=(0.8/2)^2*pi*L_engine;

D_air_inlet=0.8+2*t_fuselage;
L_air_inlet_fuselage=0.65;
V_air_inlet=2*pi*(D_air_inlet/2)^2*L_air_inlet_fuselage;
b_end=D_air_inlet/2;
a_end=b_end;

m_fuel=fuel_weight(net_thrust);
V_tot_fuel=m_fuel/800;%hypothesis: kerozen is used as the fuel.
Vf_fuel=V_tot_fuel-Vw_fuel;

%diameter of the wheels
d_nose_wheel=(5.1*(MTOW*0.12)^0.302)*10^(-2);%[m]nose wheel supports 8 to 15% of the weight.
w_nose_wheel=(0.36*(MTOW*0.12)^0.467)*10^(-2);
d_wheel=(5.1*(MTOW*(1-0.12)/2)^0.302)*10^(-2);
w_wheel=(0.36*(MTOW*(1-0.12)/2)^0.467)*10^(-2);
V_wheels=d_nose_wheel^2/4*pi*w_nose_wheel+2*d_wheel^2/4*pi*w_wheel;
V_landing_gear=V_wheels+2*(0.115/2)^2*pi*(0.8-d_wheel)+(0.115/2)^2*pi*(0.8-d_nose_wheel);

V_f_real=V_sensors+V_subsystems+V_mission+Vf_fuel+V_landing_gear+V_engine+V_air_inlet;
%%
%first estimation using statistical relations.
%MTOW and empty weight comparable to a jet trainer (resp. 3500kg and 2500kg).
A=0.333;
C=0.41;
L_first_guess=A*MTOW^C;%[m], length of the fuselage.
%%
%Volume of the fuselage, kept constant through the iterations. 
%IMPORTANT: fuselage is first considered to be cylindrical and will later
%be streamlined.
%The volume used here corresponds to the volume necessary to store all the
%payload+the subsystems+the sensors+the landing gears+fuel.
V_f=V_f_real+0.3*V_f_real;%[m^3], take a margin to ensure that the volume is still sufficient even after streamlining (30%).
f=7;%optimum fineness ratio of the fuselage to minimize drag (for a constant volume):[6;8].
%the payload is stored in a rectangle of volume V_f. The width is considered to be 1.5 larger than the height of the rectangle.
w_first_guess=sqrt(V_f/(2/3*L_first_guess));
h_first_guess=w_first_guess*2/3;
%the dimensions of the ellipse containing this rectangle:
a_first_guess=w_first_guess/2*sqrt(2);
b_first_guess=h_first_guess/2*sqrt(2);
%equivalent diameter of the ellipse:
Deq_first_guess=sqrt(4*a_first_guess*b_first_guess);
%second estimation
L_sec_guess=f*Deq_first_guess;
w_sec_guess=sqrt(V_f/(2/3*L_sec_guess));
h_sec_guess=w_sec_guess*2/3;
%the dimensions of the ellipse containing this rectangle:
a_sec_guess=w_sec_guess/2*sqrt(2);
b_sec_guess=h_sec_guess/2*sqrt(2);
%equivalent diameter of the ellipse:
Deq_sec_guess=sqrt(4*a_sec_guess*b_sec_guess);
%iterations to get the final length and diameter of the fuselage.
L_f=zeros(2,1);
Deq=zeros(2,1);
w=zeros(2,1);
h=zeros(2,1);
a=zeros(2,1);
b=zeros(2,1);
L_f(1)=L_first_guess;
L_f(2)=L_sec_guess;
w(1)=w_first_guess;
w(2)=w_sec_guess;
h(1)=h_first_guess;
h(2)=h_sec_guess;
a(1)=a_first_guess;
a(2)=a_sec_guess;
b(1)=b_first_guess;
b(2)=b_sec_guess;
Deq(1)=Deq_first_guess;
Deq(2)=Deq_sec_guess;
i=2;
while abs(L_f(i)-L_f(i-1))>0.01 && abs(Deq(i)-Deq(i-1))>0.001
    L_f(i+1)=f*Deq(i);
    w(i+1)=sqrt(V_f/(2/3*L_f(i+1)));
    h(i+1)=w(i+1)*2/3;
    a(i+1)=w(i+1)/2*sqrt(2);
    b(i+1)=h(i+1)/2*sqrt(2);
    Deq(i+1)=sqrt(4*a(i+1)*b(i+1));
    i=i+1;
end
a_val=a(i)+t_fuselage;
b_val=b(i)+t_fuselage;
L_f_val=L_f(i);
% Change the value of L_f after check of internal arrangement;

Deq_val=Deq(i);
% fprintf('The final length of the fuselage is equal to %d m \n',L_f(i));
% fprintf('The final equivalent diameter of the fuselage is equal to %d m \n',Deq(i));
fprintf('The final width of the rectangle is equal to %d m \n',w(i));
fprintf('The final height is equal to %d m \n',h(i));
% fprintf('The final dimensions of the elliptical cross-section equal to a=%d m and b=%d m\n',a(i),b(i));
%%
x_fuel=4.4;
x_engine=L_f_val-L_engine;
x_end=L_f_val;
x_air_inlet=5.2128;
b_engine=b_end;
m_a=(a_end-a_val)/(x_end-x_fuel);
p_a=a_val-m_a*x_fuel;
m_b=(2*b_engine-2*b_val)/(x_engine-x_fuel);
p_b=2*b_val-m_b*x_fuel;
syms a_x(x) b_x(x)
a_x(x)=m_a*x+p_a;
b_x(x)=m_b*x+p_b;
double(1000*a_x(x_air_inlet));
double(1000*b_x(x_air_inlet)/2);
end