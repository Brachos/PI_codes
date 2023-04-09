%Function to compute the load
function [Tx,Ty,Tz,Mx,My,Mz] = wing_load(W,Wing_loading,i,y_cg,y_ac,Mom_wing,n,Wing)
    %WING_LOAD returns the forces and moments applied at the root of the wing along each direction.
    % -----
    %
    % Syntax:
    %   [Tx,Ty,Tz,Mx,My,Mz] = wing_load(W,Wing_loading,i,Wing.aoa_fuselage,y_cg,y_ac,Mom_wing,n,Wing) returns an array with
    %   the shear stress along x, y and z [N] and the moments around x, y and z [Nm], the fuselage moment [Nm].
    %
    % Inputs:
    %   W: weight of one wing [N]
    %   Wing_loading:   Wing_loading.L: Full lift of the wing [N]    
    %                   Wing_loading.D: Full drag of the wing [N]
    %   i: the angle of attack of the wing[rad]
    %   y_ac: spanwise position of the aerodynamic center [m]
    %   y_cg: spanwise position of the center of gravity [m]
    %   Mom_wing: Wing moment [Nm].
    %   Wing:       Wing.S: Wing surface [m²]
    %               Wing.C_D: Wing drag coefficient
    %               Wing.c: Wing mean aerodynamic chord [m]
    %               Wing.C_M: Wing pitching moment coefficient [Nm]
    %               Wing.l: Wing lever arm [m]
    %               Wing.aoa: Angle of attack between the wing lever
    %               arm and the chord [rad]
    %               Wing.aoa_fuselage: Angle of attack between the wing and
    %               the fuselage center line [rad]
    %   
    % Outputs: 
    %   Tx: the shear stress along x at the wing root [N]
    %   Ty: the shear stress along y at the wing root [N]
    %   Tz: the shear stress along z at the wing root [N]
    %   Mx: the momemtum around the x-axis at the wing root [Nm]
    %   My: the momentum around the y-axis at the wing root [Nm]
    %   Mz: the momentum around the z-axis at the wing root [Nm]
    %
    %   Only half of the wing is considered.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tx = (Wing_loading.L/2)*sin(i) + 0.5*(Wing_loading.D)*cos(i);
Ty = 0;
Tz = (Wing_loading.L/2)*cos(i) + 0.5*(Wing_loading.D)*sin(i) - n*W*cos(i);
Mx = -((Wing_loading.L/2*cos(i) + Wing_loading.D/2*sin(i))*y_ac + n*W*cos(i)*y_cg);
My = Wing_loading.L/2*Wing.MAC/4 + Mom_wing/2 + (n*W*Wing.MAC/2);
Mz = (-Wing_loading.D/2*cos(i) + Wing_loading.L/2*sin(i))*y_ac - n*W*sin(i)*y_cg;
end