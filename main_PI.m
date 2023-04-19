clear
MTOW = 4082; %first guess
cg = 4.444; %first guess
drag = 3278; %first guess
error = 1;
prop_lift = 0.844;
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