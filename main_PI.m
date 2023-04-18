clear
MTOW = 3050; %first guess
cg = 3.5; %first guess
error = 1;
prop_lift = 0.9;
while error >  1e-4
    [new_MTOW, new_cg, new_prop_lift WING, V_TAIL, RUDDER] = Aircraft(MTOW, cg,prop_lift);
    error = max([abs(new_MTOW - MTOW)/MTOW abs(new_cg(1) - cg)/cg abs(new_prop_lift-prop_lift)/prop_lift]);
    MTOW = new_MTOW;
    cg = new_cg(1);
    prop_lift = new_prop_lift;
end