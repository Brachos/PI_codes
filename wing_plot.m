%% Function that plot all graphs and figures relative to the wing design

function [] = wing_plot(WING,FUSELAGE)

    % Color
        gray = [.8 .8 .8];
    
    
    % Wing geometry (including flaps and ailerons)
        % Draw
        c_root = WING.cw_root;
        c_tip = WING.cw_tip;
        b = WING.bw;
        sweep = WING.sweep;
        R_fus = FUSELAGE.a_el;
        L_HL = atan2(0.45*(c_tip-c_root)+b/2*tan(sweep),b/2); % Sweep (HL)
        dx_tip = (c_root-c_tip)/4 + b*tan(sweep)/2;
        wing_x = [c_root c_root-dx_tip c_root-dx_tip-c_tip 0]; 
        wing_x(5:8) = wing_x;
        wing_y = [0 -b/2 -b/2 0 0 b/2 b/2 0];
        SL = [0.75*c_root c_root-dx_tip-0.25*c_tip;0 -b/2]; % Sweep line
        SHL = [0.75*c_root 0.75*c_root;0 -b/2]; % Sweep horizon line
        HL = [0.3*c_root 0.3*c_root-tan(L_HL)*b/2;0 -b/2]; % Hinge line
        HHL = [0.3*c_root 0.3*c_root;0 -b/2]; % Hinge horizon line
        pente_HL = (HL(1,1)-HL(1,2))/(HL(2,1)-HL(2,2));
        pente_tr = (0-wing_x(3))/(b/2);
        fus_l = [-R_fus -R_fus; -0.3*c_root 1.3*c_root];
        fus_r = [R_fus R_fus; -0.3*c_root 1.3*c_root];
        flap_yl = [-R_fus-0.1 -R_fus-0.1-b/2*0.4 -R_fus-0.1-b/2*0.4 -R_fus-0.1 -R_fus-0.1];
        flap_yr = -flap_yl;
        fx1 = 0.3*c_root+flap_yl(1)*pente_HL; % First point of the flap
        fx2 = 0.3*c_root+flap_yl(2)*pente_HL;
        fx3 = flap_yl(2)*pente_tr;
        fx4 = flap_yl(1)*pente_tr;
        flap_x = [fx1 fx2 fx3 fx4 fx1];
        aileron_yl = [-R_fus-0.1-b/2*0.4 -b/2+0.1 -b/2+0.1 -R_fus-0.1-b/2*0.4 -R_fus-0.1-b/2*0.4];
        aileron_yr = - aileron_yl;
        ax1 = 0.3*c_root+aileron_yl(1)*pente_HL; % First point of the flap
        ax2 = 0.3*c_root+aileron_yl(2)*pente_HL;
        ax3 = aileron_yl(2)*pente_tr;
        ax4 = aileron_yl(1)*pente_tr;
        aileron_x = [ax1 ax2 ax3 ax4 ax1];
        % Mesures
        m_b = [-b/2 b/2;-c_root -c_root];
       
    
    figure
    plot(wing_y,wing_x,'k','linewidth',0.7)
    hold on
    plot(SL(2,:),SL(1,:),'--k','linewidth',0.5)
    plot(SHL(2,:),SHL(1,:),'--k','linewidth',0.5)
    plot(HL(2,:),HL(1,:),'--k','linewidth',0.5)
    plot(HHL(2,:),HHL(1,:),'--k','linewidth',0.5)
    plot(fus_l(1,:),fus_l(2,:),'--','Color',gray,'linewidth',0.3)
    plot(fus_r(1,:),fus_r(2,:),'--','Color',gray,'linewidth',0.3)
    plot(flap_yl,flap_x,'k','linewidth',0.7)
    plot(flap_yr,flap_x,'k','linewidth',0.7)
    plot(aileron_yl,aileron_x,'k','linewidth',0.7)
    plot(aileron_yr,aileron_x,'k','linewidth',0.7)
    axis equal
    axis off
    
    % Drag polar
    
    
    
    % CL vs AOA