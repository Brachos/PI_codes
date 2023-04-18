clear
MTOW = 3050; %first guess
cg = 3.5; %first guess
error = 1;
while error >  1e-4
    [new_MTOW, new_cg, WING, V_TAIL, RUDDER] = Aircraft(MTOW, cg);
    error = max(abs(new_MTOW - MTOW)/MTOW, abs(new_cg(1) - cg)/cg);
    MTOW = new_MTOW;
    cg = new_cg(1);
end