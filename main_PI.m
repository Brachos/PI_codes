clear
MTOW = 3213.2; %first guess
cg = 4.3373; %first guess
drag = 2199; %first guess
error = 1;
prop_lift = 0.8004;

while error >  1e-4

    [new_MTOW, new_cg, new_prop_lift, WING, V_TAIL, RUDDER, FUSELAGE, WEIGHT, PARAM,COST, LongDeriv, LatDeriv, Long_modes, LatModes] = Aircraft(MTOW, cg,prop_lift, drag);
    [DRAG] = Drag(WING, V_TAIL, FUSELAGE, PARAM);
    new_drag = 1/2 * DRAG.CDTotal * PARAM.rho * WING.Sw * PARAM.V_c^2;
    error = max([abs(new_MTOW - MTOW)/MTOW abs(new_cg(1) - cg)/cg abs(new_prop_lift-prop_lift)/prop_lift abs(new_drag - drag)/drag]);
    MTOW = new_MTOW;
    cg = real(new_cg(1));
    prop_lift = new_prop_lift;
    drag = new_drag;
end
%
g = 9.81;
lift_to_drag=MTOW*g/drag;

% Verifications
if drag < 3500
    D = 'OK';
else 
    D = 'NOT OK';
end
fprintf('Drag is %s\n',D);


if lift_to_drag > 10
    LTD = 'OK';
else 
    LTD = 'NOT OK';
end
fprintf('Lift to drag ratio is %s\n',LTD);

if WING.bw < 10
    span = 'OK';
else 
    span = 'NOT OK';
end
fprintf('Wing span is %s\n',span);

if WING.Sw < 10
    Area = 'OK';
else 
    Area = 'NOT OK';
end
fprintf('Wing surface is %s\n',Area);
%coucou
if prop_lift > 0.8 && prop_lift < 1
    PL = 'OK';
else 
    PL = 'NOT OK';
end

fprintf('Prop_lift is %s\n',PL);

if PARAM.kf < 0.2 && PARAM.kf > 0.05
    KF = 'OK';
else 
    KF = 'NOT OK';
end
fprintf('Static margin for aircraft with fuel but no payload is %s\n',KF);


if PARAM.k < 0.2 && PARAM.k > 0.05
    K = 'OK';
else 
    K = 'NOT OK';
end
fprintf('Static margin for aircraft fully loaded is %s\n',K);

if PARAM.k2 < 0.2 && PARAM.k2 > 0.05
    K2 = 'OK';
else 
    K2 = 'NOT OK';
end
fprintf('Static margin for aircraft empty is %s\n',K2);

if PARAM.kp < 0.2 && PARAM.kp > 0.05
    KP = 'OK';
else 
    KP = 'NOT OK';
end
fprintf('Static margin for aircraft with payload but without fuel is %s\n',KP);


%% Break-even analysis 
% [new_MTOW, new_cg, new_prop_lift, WING, V_TAIL, RUDDER, FUSELAGE, WEIGHT, PARAM,COST] = Aircraft(MTOW, cg,prop_lift, drag,N_old);
% Total_fixed_cost = COST.Ccert;
% Unit_variable_cost = Total_fixed_cost/N_old;
% N_vector = 1:1:1000;
% Variable_cost = [];
% unite_sales_price = [500000 600000 700000 800000 900000 1000000];
% saling_price_per_production = [];
% for j=1:length(unite_sales_price)
%     for i=1:length(N_vector)
%         saling_price_per_production(j,i) = N_vector(i)*unite_sales_price(j);
%     end
% end
% for i =1:length(unite_sales_price)
%     NBE(i) = Total_fixed_cost/(unite_sales_price(i) - Unit_variable_cost);
% end
% for i=1:length(N_vector)
%     [new_MTOW, new_cg, new_prop_lift, WING, V_TAIL, RUDDER, FUSELAGE, WEIGHT, PARAM,COST] = Aircraft(MTOW, cg,prop_lift, drag,N_old);
%     Variable_cost(i) = COST.CMFG + COST.Cqc + COST.Cmat + COST.Cpp + 15000;
%     Fixed_and_variable(i) = Total_fixed_cost + Variable_cost(i);
% end
% show = 1;
% if show == 1
%     figure;
%     a = plot(N_vector,Fixed_and_variable,'LineWidth',2)
%     hold on 
%     plot(N_vector,saling_price_per_production(1,:),'LineWidth',2)
%     plot(N_vector,saling_price_per_production(2,:),'LineWidth',2)
%     plot(N_vector,saling_price_per_production(3,:),'LineWidth',2)
%     plot(N_vector,saling_price_per_production(4,:),'LineWidth',2)
%     plot(N_vector,saling_price_per_production(5,:),'LineWidth',2)
%     plot(N_vector,saling_price_per_production(6,:),'LineWidth',2)
% end
