rho = 0.46;
V = speed(30000,0.7);
L = 8.07;
L = 1.2;
mu = 1.453e-5;
Re = rho*V*L/mu;
AOA = 0*pi/180;

% Dimensions
c_root_in  = 8.07;
c_root_out = c_root_in - 1.08;
c_tip = 1.5;
b1 = 0.68;
b2 = 0.74;
tap1 = c_root_out/c_root_in;
tap1 = c_tip/c_root_out;
S1 = (c_root_in + c_root_out)*0.34;
S2 = (c_root_out + c_tip)*0.37;
S = S1 + S2;
AR1 = b1^2/S1;
AR2 = b2^2/S2;

% Lift Coefficient of the 2 parts of the fuselage
beta = sqrt(1-0.7^2);
Cte = (c_root_out-c_root_in)/4 + 1.08;
sweep1 = atan2(Cte,0.34);
Cte = (c_tip-c_root_out)/4 + 1.18;
sweep2 = atan2(Cte,0.37);

L_beta_1 = atan(tan(sweep1)/beta); % Angle to graphically find x_AC
L_beta_2 = atan(tan(sweep2)/beta);

% CL and CD
cl_alpha  = (0.5384-0.1873)*180/(2*pi); % [1/rad]
alpha_l0  = -2.04*pi/180;
alpha_01  = -0.23;     % See L5 - P30
theta_tip = 0*pi/180; % [rad] Twist angle
alpha_L0  = alpha_l0 + alpha_01*theta_tip;
k  = beta*cl_alpha/(2*pi);

a_1 = 2*pi/(2/(beta*AR1)+sqrt((1/(k*cos(L_beta_1)))^2+(2/(beta*AR1))^2))/beta;
a_2 = 2*pi/(2/(beta*AR2)+sqrt((1/(k*cos(L_beta_2)))^2+(2/(beta*AR2))^2))/beta;
CL_alpha_1 = a_1;
CL_alpha_2 = a_2;

CL_1 = a_1*(AOA-alpha_L0);
CL_2 = a_2*(AOA-alpha_L0);

CL = (CL_1*S1+CL_2*S2)/S;
L  = 0.5*rho*S*CL*V^2;
m = L/9.81;

% MAC
syms y
c    = (1-2*y/1.42)*c_root_in + (2*y/1.42)*c_tip;
c_AC = double(2/6.8*int(c^2,0,1.42/2));