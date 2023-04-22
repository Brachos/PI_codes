function [CL_max_TO, CL_max_L, fb_ratio, S_flap] = flaps(WING,FUSELAGE)

    % fb_ratio = flap span ratio
    % S_flap = surface of flap
    
    c_root = WING.cw_root;
    c_tip = WING.cw_tip;
    b = WING.bw;
    S = WING.Sw;
    sweep = WING.sweep;
    R_fus = FUSELAGE.a_el;
    dx_tip = (c_root-c_tip)/4 + b*tan(sweep)/2;
    L_HL = atan2(0.45*(c_tip-c_root)+b/2*tan(sweep),b/2); % Sweep (HL)
    
    pente_lead  = dx_tip/(b/2);
    pente_trail = -(c_root-dx_tip-c_tip)/(b/2);
    y_flap = [-R_fus-0.1 -R_fus/2-b/4];
    xf1 = c_root + pente_lead*y_flap(1);
    xf2 = pente_trail*y_flap(1);
    xf3 = c_root + pente_lead*y_flap(2);
    xf4 = pente_trail*y_flap(2);
    
    S_flap = (xf1-xf2+xf3-xf4)*(y_flap(1)-y_flap(2));
    fb_ratio = (y_flap(1)-y_flap(2))/(b/2);
    
    CL_max_base = 0.95*cos(sweep)*(1.75+1.55)/2*1.1; % Graph L5-S32
    DCL_max_TO  = 0.9*cos(L_HL)*(S_flap/S)*0.95; % Graph L5-S33 with single slotted, Angle = 20 deg 
    DCL_max_L   = 0.9*cos(L_HL)*(S_flap/S)*1.5;  % Angle = 40 deg 
    CL_max_TO = CL_max_base + DCL_max_TO;
    CL_max_L  = CL_max_base + DCL_max_L;
end