clear
MTOW = 4082; %first guess
cg = 4.444; %first guess
drag = 3278; %first guess
error = 1;
prop_lift = 0.7612;
while error >  1e-4
    [new_MTOW, new_cg, new_prop_lift, WING, V_TAIL, RUDDER, FUSELAGE, WEIGHT, PARAM] = Aircraft(MTOW, cg,prop_lift, drag);
    [DRAG] = Drag(WING, V_TAIL, FUSELAGE, PARAM);
    new_drag = 1/2 * DRAG.CDTotal * PARAM.rho * WING.Sw * PARAM.V_c^2;
    error = max([abs(new_MTOW - MTOW)/MTOW abs(new_cg(1) - cg)/cg abs(new_prop_lift-prop_lift)/prop_lift abs(new_drag - drag)/drag]);
    MTOW = new_MTOW;
    cg = real(new_cg(1));
    prop_lift = new_prop_lift;
    drag = new_drag;
end
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

if prop_lift > 0.85 && prop_lift < 1
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
