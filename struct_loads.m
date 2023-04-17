function[Tx,Ty,Tz,Mx,My,Mz] = struct_loads(Fuselage,Wing,Tail,Engine,Sensors,Rear_land_gear,Payload,Empennage,Fin,aoa,n, Lift_empennage, F_fin, M_fus)
    %STRUCT_LOADS  returns the forces and moments applied at the rear of the fuselage along each direction.
    %The structural loads are computed with respect to the back of the wing
    %and with respect to an intermediate section which are called section A.
    % -----
    %
    % Syntax:
    %   [Tx,Ty,Tz,Mx,My,Mz] = struct_loads(Fuselage,Wing,Tail,Engine,Sensors,Rear_land_gear,Payload,Empennage,aoa,n, Lift_empennage, F_fin, M_fus)
    %   returns an array with the shear stress along x, y and z [N] and the moments around x, y and z [Nm], the fuselage moment [Nm] computed at 
    %   the rear of the fuselage.
    %
    % Inputs:
    %
    % Fuselage:     Fuselage.x_min: x-position of the minimum cross-section (at tail)[m]
    %               Fuselage.a_min: Minimum semi major axis of the elliptical fuselage (at the tail) [m]
    %               Fuselage.b_min: Minimum semi minor axis of the elliptical fuselage (at the tail) [m]
    %               Fuselage.x_min: x-position of the maximum cross-section (at tail)[m]
    %               Fuselage.a_max: Semi major axis of the fuselage at section A [m]
    %               Fuselage.b_max: Semi minor axis of the fuselage at section A [m]
    %               Fuselage.x_cs: x-position of the considered cross-section [m]
    %               Fuselage.a: Semi major axis of the considered cross-section [m]
    %               Fuselage.b: Semi-minor axis of the considered cross-section [m]
    %               Fuselage.A_h: Area of the cross-sections[m²]
    %               Fuselage.L: Total length of the fuselage [m]
    %               Fuselage.W: Total weight of the fuselage [N]
    %   
    %   Wing:       Wing.S: Wing surface [m²]
    %               Wing.W: Wing weight [N]
    %               Wing.AR: Wing aspect ratio [-]
    %               Wing.MAC: Wing mean aerodynamic chord [m]
    %               Wing.l: Wing lever arm [m]
    %               Wing.aoa: Angle of attack between the wing lever
    %               arm and the chord [rad]
    %               Wing.LE: Position of the leading edge [m]
    %               Wing.root_chord: Wiing root chord [m]
    %               Wing.aoa_fuselage: Angle of attack between the wing and
    %               the fuselage center line [rad]
    %
    %   Tail:       Tail.cg: Position of the cg of the tail [m]
    %               Tail.W: Weight of the v-tail [N]
    %
    %   Engine:     Engine.W: Engine weight [N]
    %               Engine.cg: Center of gravity of the engine [m]
    %
    %   Sensors:    Sensors.W: Sensors weight [N]
    %               Sensors.cg: Center of gravity of the sensors [m]
    %
    %   Rear landing gear: Rear_land_gear.W: Rear landing gear weight [N]
    %               Rear_land_gear.cg: Center of gravity of the Rear landing gear [m]
    %
    %   Payload:    Payload.W: Payload weight [N]
    %               Payload.cg: Center of gravity of the Payload [m]
    %
    %   Empennage:  Empennage.S: Empennage surface [m²]
    %               Empennage.W: Empennage weight [N]
    %               Empennage.MAC: Empennage mean aerodynamic chord [m]
    %               Empennage.b: Empennage span [m]
    %               Empennage.l: Empennage lever arm [m]
    %               Empennage.T: Thrust placed on Empennage [N]
    %               Empennage.l_T: Empennage thrust lever arm [m]
    %               Empennage.aoa: Angle of attack between the empennage lever
    %               arm and the chord [rad]
    %               Empennage.aoa_T: Angle of attack between the thrust lift lever
    %               arm and the chord [rad]
    %               Empennage.ac: Aerodynamic centre of the empennage [m]
    %               Empennage.aoa_fuselage: Angle of attack between the empennage and
    %               the fuselage center line [rad]
    %
    %   Fin:        Fin.AR: Fin aspect ratio [-]
    %               Fin.S: Fin surface [m²]
    %               Fin.l: Fin lever arm [m]
    %
    %   aoa: angle of attack [rad]
    %   n: load [g]
    %   Lift_empennage: Empennage load [N]
    %   F_fin: Fin load [N]
    %   M_fus: Fuselage moment [Nm]
    %
    % Outputs:
    %   Tx: the shear stress along x at the rear of the fuselage [N]
    %   Ty: the shear stress along y at the rear of the fuselage [N]
    %   Tz: the shear stress along z at the rear of the fuselage [N]
    %   Mx: the momemtum around the x-axis at the rear of the fuselage [Nm]
    %   My: the momentum around the y-axis at the rear of the fuselage [Nm]
    %   Mz: the momentum around the z-axis at the rear of the fuselage [Nm]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_A = [Wing.LE+Wing.root_chord; Wing.LE+Wing.root_chord-(Fuselage.x_cs(2)-Fuselage.x_cs(1))];%x_A is the distance between the nose of the aircraft and section A.
L_rear = [Fuselage.L-x_A(1); Fuselage.L-x_A(2)];%length of the fuselage after the back of the wing
Tail_cg_A = [Tail.cg-x_A(1); Tail.cg-x_A(2)];   %position of the center of gravity of the
                                                %tail with respect to the back of the wings

%Computation of the skin area of the fuselage, considering it is uniformly
%tapered, after the wing
f1 = @(t1) sqrt(Fuselage.a_max^2 * (cos(t1)).^2 + Fuselage.b_max^2 * (sin(t1)).^2);
perim_max =  integral(f1,0,2*pi);
f2 = @(t2) sqrt(Fuselage.a_min^2 * (cos(t2)).^2 + Fuselage.b_min^2 * (sin(t2)).^2);
perim_min =  integral(f2,0,2*pi);
f3 = @(t3) sqrt(Fuselage.a(2)^2 * (cos(t3)).^2 + Fuselage.b(2)^2 * (sin(t3)).^2);
perim_inter =  integral(f3,0,2*pi);
area_skin = [1/2*(perim_min+perim_max)*L_rear(1); 1/2*(perim_min+perim_inter)*L_rear(2)];

%self-weight of the fuselage per meter 
W_rear = Fuselage.W.*(L_rear./Fuselage.L);  %weight of the empty fuselage after the wing
weight_per_meter_min = W_rear.* perim_min./area_skin;
weight_per_meter_max = [W_rear(1)*perim_max/area_skin(1); W_rear(2)*perim_inter/area_skin(2)];

%the self-weight induces some shear forces (SF) and bending moments (BM)
%SF and BM are calculated by equilibrium ('MNT'), it gives here a triangle
%repartition q1 and a constant linear force q2 
q1 = [linspace(weight_per_meter_max(1)- weight_per_meter_min(1),0, 100); linspace(weight_per_meter_max(2)- weight_per_meter_min(2),0, 100)];
q2 = weight_per_meter_min;

%Resultants of the 2 distribution of forces
Q1 = [L_rear(1)*q1(1, 1)/2; L_rear(2)*q(2, 1)/2];
Q2 = L_rear.*q2;

SF = [n*(Q1(1)+Q2(1)+Tail.W+Engine.W+Sensors.W+Rear_land_gear.W+Payload.W);
    n*(Q1(2)+Q2(2)+Tail.W+Engine.W+Sensors.W+Rear_land_gear.W+Payload.W)];
BM = [n*cos(aoa-Wing.aoa_fuselage)*(Tail_cg_A(1)*Tail.W + 1/3*L_rear(1)*Q1(1) + 1/2*L_rear(1)*Q2(1) + Engine.W*(Engine.cg-x_A(1)) + Sensors.W*(Sensors.cg-x_A(1)) + Rear_land_gear.W*(Rear_land_gear.cg-x_A(1)) + Payload.W*(Payload.cg-x_A(1)));
    n*cos(aoa-Wing.aoa_fuselage)*(Tail_cg_A(2)*Tail.W + 1/3*L_rear(2)*Q1(2) + 1/2*L_rear(2)*Q2(2) + Engine.W*(Engine.cg-x_A(2)) + Sensors.W*(Sensors.cg-x_A(2)) + Rear_land_gear.W*(Rear_land_gear.cg-x_A(2)) + Payload.W*(Payload.cg-x_A(2)))];

%Resultant forces in the section after the wing
Tx = [(-SF(1)+Lift_empennage)*sin(aoa-Wing.aoa_fuselage); (-SF(2)+Lift_empennage)*sin(aoa-Wing.aoa_fuselage)];
Ty = -F_fin.*[1; 1];
Tz = [(SF(1)-Lift_empennage)*cos(aoa-Wing.aoa_fuselage);(SF(2)-Lift_empennage)*cos(aoa-Wing.aoa_fuselage)];
Mx = -M_fus.*[1; 1];   %it takes into account the fin and the empennage
My = [BM(1)-Lift_empennage*(Empennage.ac-x_A(1))*cos(aoa-Wing.aoa_fuselage); BM(2)-Lift_empennage*(Empennage.ac-x_A(2))*cos(aoa-Wing.aoa_fuselage)];
Mz = [F_fin*(Fin.ac-x_A(1)); F_fin*(Fin.ac-x_A(2))];

end