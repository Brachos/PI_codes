function [Deq_val,a_val,b_val,L_f_val,V_f]=fuselage_design(MTOW,Vw_fuel)
%volume of the payload
V_sensors=1;%[m^3]
V_subsystems=1;%[m^3]
V_mission=0.15;%[m^3]
m_fuel=2221;%[kg]
V_tot_fuel=m_fuel/800;%hypothesis: kerozen is used as the fuel.
Vf_fuel=V_tot_fuel-Vw_fuel;
%diameter of the wheels
d_nose_wheel=(5.1*(MTOW*0.12)^0.302)*10^(-2);%[m]nose wheel supports 8 to 15% of the weight.
w_nose_wheel=(0.36*(MTOW*0.12)^0.467)*10^(-2);
d_wheel=(5.1*(MTOW*(1-0.12)/2)^0.302)*10^(-2);
w_wheel=(0.36*(MTOW*(1-0.12)/2)^0.467)*10^(-2);
V_wheels=d_nose_wheel^2/4*pi*w_nose_wheel+2*d_wheel^2/4*pi*w_wheel;
V_landing_gear=V_wheels;
V_f_real=V_sensors+V_subsystems+V_mission+Vf_fuel+V_landing_gear;
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
V_f=V_f_real+0.3*V_f_real;%[m^3], take a margin to ensure that the volume is still sufficient even after streamlining.
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
a_val=a(i);
b_val=b(i);
L_f_val=L_f(i);
Deq_val=Deq(i);
% fprintf('The final length of the fuselage is equal to %d m \n',L_f(i));
% fprintf('The final equivalent diameter of the fuselage is equal to %d m \n',Deq(i));
 fprintf('The final width of the rectangle is equal to %d m \n',w(i));
 fprintf('The final height is equal to %d m \n',h(i));
% fprintf('The final dimensions of the elliptical cross-section equal to a=%d m and b=%d m\n',a(i),b(i));
end