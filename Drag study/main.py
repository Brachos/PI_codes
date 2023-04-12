import PlaneWings, tailDrag, engine, fuselageLiftDrag, InterferenceEffectAdCorrections, \
    Protuberances_surfaceImperfections
import math
import computing

""" ------------Aircraft parameters------------"""

# wing parameters:
S =  # Wing surface area
S_net =  # surface with fuselage substracted, A-3.1
mach =  # flight Mach number?

# value found in "C_values.png", there is a function in "computing.py" to get the abciss of this graph
C_1 =
C_2 =
C_3 =

taper_ratio =
sweepback_angle =  # sweepback angle of quarter chord
C_L =  # 3D lift coefficient
C_L_i =  # design lift coefficient (mean of the 2D coef on the whole wing
C_L_max =  # E-4.4
t =  # mean thickness
c =  # mean chord length
shape =  # shape factor of wing, wing_shape_factor.png
C_F =  # skin friction wing, found in "skinFriction.png", function in"computing.py" to help find it
Lambda_half =  # sweepback angle, of half chord?
aspect_ratio =  # of wing
epsilon_t =  # angle of twist
C_o_l =  # factor found in figure "wingTwistFactors.png"
C_l_l =  # actor found in figure "wingTwistFactors.png"

# engine parameters:
skin_friction_engine =  # found in "skinFriction.png", function in"computing.py" to help find it
S_n_wet_engine =  # wet surface of the engine

l_p =  # plug length
D_p =  # plug diameter

M_g =  # engine type specification or H-20

# pylon properties, not used
"""py_thickness =  # relative thickness of pylon in streamwise direction
py_c =  # ?
Lambda_py =  #
S_py_wet =  # wetted area of pylon"""

# fuselage parameters
C_F_fuse =  # skin friction fuselage, found in "skinFriction.png", function in"computing.py" to help find it
S_f_wet_fuse =  # wetted area of fuselage
A_c =  # cross section of fuselage (can take the mean i guess)
l_n =  # fuselageNomenclature.png
l_c =  # fuselageNomenclature.png
C_D_beta =  # C_D_beta.png
C_L_0 =  # C_L at AOA of 0 relative to fuselage
C_L_alpha =  # lift curve slope of aircraft (rad^-1)
V_f =  # volume of the fuselage

# tail parameters:
C_F_tail =  # found in "skinFriction.png", function in"computing.py" to help find it
t_tail =  # average thickness of tail
c_tail =  # average chord of tail
Lambda_half_tail =  # sweep angle at half chord
S_tail =  # gross area of tail
C_L_h =  # tailplane lift, E-44
A_h =  # aspect ratio of tail

# imperfection parameters:
flaps_dragArea =  # flaps.png
imperfection_dragArea =  # imperfections.png

# interferences and corrections parameters:


t_r =  # thickness of wing at root
C_ci =  # total circumferential Cl.length for both wing halves of the wing fuselage intersection line at which the boundary layers interact, 4.5 * the root chord for a midwing
c_r =  # chord at root
D_F =  # Fuselage diameter
beta =  # rear fuselage upsweep angle (can be negative?)

d_2 =  # F-53, D_parameters.png
deps_dC_L =  # E-10.1

# computed parameters:
compressibility_factor = math.sqrt(1 - mach ** 2)
corrected_sweepback_angle = math.atan(math.tan(sweepback_angle) / compressibility_Factor)
S_p_wet = computing.plug_wetted_area(l_p, D_p)  # B-4, plug wetted area

"""----------------------------Drag calculation----------------------------"""

# wing
untwistWing = PlaneWings.untwistedPlaneWings(C_1, C_2, C_3, taper_ratio, corrected_sweepback_angle, C_L,
                                             aspect_ratio)
wing = PlaneWings.wing3D(C_L_max, t, c, C_L, C_L_i, C_F, shape, Lambda_half, S,
                         S_net)  # typically between .006 and .01

wingTwist = PlaneWings.wingTwist(epsilon_t, C_o_l, C_L, C_l_l)

wingDrag = untwistWing + wing + wingTwist

# Engine
fanCowlingDrag = engine.FanCowlingDragArea_simple(skin_friction_engine, S_n_wet_engine)
gasGenCowlingDrag = engine.gasGeneratorCowlingDragComponent_Estimation(skin_friction_engine, S_p_wet, M_g, mach)
# pyDrag = pylonDrag(skin_friction_engine, py_thickness, py_c, Lambda_py, S_py_wet) #we don't use pylon, as the engine is attached to the body

engineDragArea = fanCowlingDrag + gasGenCowlingDrag + pyDrag
engineDrag = engineDragArea / S

# tail
tailBasic = tailDrag.BasicDrag(C_F_tail, t_tail, c_tail, Lambda_half_tail, S_tail)
tailLiftDrag = tailDrag.DragDueToLift(C_L_h, Lambda_h_tail, A_h, S_h_tail)

tailDragArea = tailBasic + tailLiftDrag
tailDrag = tailDragArea / S

# FuselageDrag

fuselageLiftDrag = fuselageLiftDrag.fuselag(C_F_fuse, S_f_wet_fuse, A_c, l_n, l_c, C_D_beta)
vortexFuseDrag = fuselageLiftDrag.VortexDragByFuselageLift(C_L, C_L_0, C_L, C_L_alpha, V_f)
fuselageDragArea = fuselageLiftDrag + vortexFuseDrag
fuselageDrag = fuselageDragArea / S

# interferance and corrections

viscousInter = InterferenceEffectAdCorrections.visoucsInterference(C_F_fuse, t_r, C_ci)
wingCorr = InterferenceEffectAdCorrections.wingBaseCorrection(C_F, c_r, D_F)
upWash = InterferenceEffectAdCorrections.upwashCorrection(beta, sweepback_angle, A, d_2, C_L, C_L_alpha)
wingTailInter = InterferenceEffectAdCorrections.wingTailplaneInterference_vortexInduced(C_L_h, C_L, deps_dC_L,
                                                                                        aspect_ratio,
                                                                                        S_tail)

InterfDragArea = viscousInter + wingCorr + upWash + wingTailInter
InterfDrag = InterfDragArea / S

# Protuberances and surface imperfections
protrudingIntake = Protuberances_surfaceImperfections.protrudingIntake(Area)
imperfectionDragArea = protrudingIntake + flaps_dragArea + imperfection_dragArea
imperfectionDrag = imperfectionDragArea / S

drag = wingDrag + engineDrag + tailDrag + fuselageDrag + InterfDrag + imperfectionDrag
print(Drag)
