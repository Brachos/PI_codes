import PlaneWings, tailDrag, engine, fuselageLiftDrag, InterferenceEffectAdCorrections, \
    Protuberances_surfaceImperfections
import math
import computing

""" ------------Aircraft parameters------------"""

# wing parameters:
S = 12.1  # Wing surface area
S_net = 8.06  # surface with fuselage substracted, A-3.1
mach = 0.7  # flight Mach number

# value found in "C_values.png", there is a function in "computing.py" to get the abciss of this graph
C_1 = .5
C_2 = .12
C_3 = .38

taper_ratio = 0.3
sweepback_angle = 36.6 / 180 * math.pi  # sweepback angle of quarter chord
C_L = 0.7  # 3D lift coefficient
C_L_i = 0.61  # design lift coefficient (mean of the 2D coef on the whole wing)
C_L_max = 1.665  # E-4.4
t = .157  # mean thickness
c = 1.126  # mean chord length
shape = .256  # shape factor of wing, wing_shape_factor.png
C_F = .003  # skin friction wing, found in "skinFriction.png", function in"computing.py" to help find it
Lambda_half = 36.6 / 180 * math.pi  # sweepback angle, of half chord
aspect_ratio = 7  # of wing
epsilon_t = 0  # angle of twist
# C_o_l = 1  # factor found in figure "wingTwistFactors.png"
# C_l_l = 1  # actor found in figure "wingTwistFactors.png"

# engine parameters:
skin_friction_engine = 0.05  # found in "skinFriction.png", function in"computing.py" to help find it
S_n_wet_engine = 1.23  # wet surface of the engine


# fuselage parameters
C_F_fuse = .0025  # skin friction fuselage, found in "skinFriction.png", function in"computing.py" to help find it
S_f_wet_fuse = 81.62  # wetted area of fuselage
A_c = 1.2  # cross-section of fuselage
l_n = 2.24  # fuselageNomenclature.png
l_c = 1.3  # fuselageNomenclature.png
C_D_beta = .1  # C_D_beta.png
C_L_alpha = 5.45  # lift curve slope of aircraft (rad^-1)
C_L_0 = C_L_alpha * 0.043  # C_L at AOA of 0 relative to fuselage
# c_l_alpha = 1  # lift curve slope of wing 2D
V_f = 6.54  # volume of the fuselage

# tail parameters:
C_F_tail = .0023  # found in "skinFriction.png", function in"computing.py" to help find it
t_tail = 0.1  # average thickness of tail
c_tail = 0.8  # average chord of tail
Lambda_half_tail = 1  # sweep angle at half chord
S_tail = 2.76  # gross area of tail
C_L_h = 0.11  # tailplane lift, E-44
A_h = 3  # aspect ratio of tail

# imperfection parameters:

Area = 0.14325  # air intake area

# interferences and corrections parameters:
t_r = 0.2  # thickness of wing at root
C_ci = 4.5 * 1.68  # total circumferential length for both wing halves of the wing fuselage intersection line at which the boundary layers interact, 4.5 * the root chord for a midwing
c_r = 1.68  # chord at root
D_F = 1.2  # Fuselage diameter
beta = 0  # rear fuselage upsweep angle

# d_2 = 1  # F-53, D_parameters.png
r = 0.73  # longitudinal position of horizontal stabilizer fig E-13
m = 0.05  # vertical position of the horzontal stabilizer fig E-13

# computed parameters:
compressibility_factor = math.sqrt(1 - mach ** 2)
corrected_sweepback_angle = math.atan(math.tan(sweepback_angle) / compressibility_factor)
# S_p_wet = computing.plug_wetted_area(l_p, D_p)  # B-4, plug wetted area
E = 1 + 2 * taper_ratio / (aspect_ratio * (1 + taper_ratio))
C_L_W_alpha = 0.085 / 180 * math.pi  # .995 * c_l_alpha / (E + c_l_alpha / (math.pi * aspect_ratio))  # lift curve slope, E-3 and E-4 #E-5
deps_dC_L = C_L_W_alpha / (math.pi * aspect_ratio * (taper_ratio * r) ** (.25) * (1 + math.fabs(m)))  # E-10.1
"""----------------------------Drag calculation----------------------------"""


# wing
def wing_Drag():
    untwistWing = PlaneWings.untwistedPlaneWings(C_1, C_2, C_3, taper_ratio, corrected_sweepback_angle, C_L,
                                                 aspect_ratio)
    wing = PlaneWings.wing3D(C_L_max, t, c, C_L, C_L_i, C_F, shape, Lambda_half, S,
                             S_net)  # typically between .006 and .01

    # wingTwist = PlaneWings.wingTwist(epsilon_t, C_o_l, C_L, C_l_l)

    wingDrag = untwistWing + wing  # + wingTwist
    return wingDrag


# Engine
def engine_Drag():
    fanCowlingDrag = engine.FanCowlingDragArea_simple(skin_friction_engine, S_n_wet_engine)
    # gasGenCowlingDrag = engine.plugDragArea(skin_friction_engine, S_p_wet, M_g, mach) #no plug
    # pyDrag = pylonDrag(skin_friction_engine, py_thickness, py_c, Lambda_py, S_py_wet) #we don't use pylon, as the engine is attached to the body

    engineDragArea = fanCowlingDrag  # + gasGenCowlingDrag  # + pyDrag
    engineDrag = engineDragArea / S
    return engineDrag


# tail
def tail_Drag():
    tailBasic = tailDrag.BasicDrag(C_F_tail, t_tail, c_tail, Lambda_half_tail,
                                   S_tail)  # To account for interference with theairplane, Sh rnay be assurned as the grosshorizontal tailplane area, including parts covered by the fuselage or vertical tailplane (cf. Section F-4.4 . )
    tailLiftDrag = tailDrag.DragDueToLift(C_L_h, Lambda_half_tail, A_h, S_tail)
    tailVortex = tailDrag.tailVortexDrag(C_L_h, S_tail, A_h)
    tailDragArea = tailBasic + tailLiftDrag + tailVortex
    tail_Drag = tailDragArea / S
    return tail_Drag


# FuselageDrag
def fuselage_Drag():
    fuselage_LiftDrag = fuselageLiftDrag.fuselag(C_F_fuse, S_f_wet_fuse, A_c, l_n, l_c, C_D_beta)
    print(fuselage_LiftDrag)
    vortexFuseDrag = fuselageLiftDrag.VortexDragByFuselageLift(C_L, C_L_0, C_L_alpha, V_f)
    print(vortexFuseDrag)
    fuselageDragArea = fuselage_LiftDrag + vortexFuseDrag
    fuselageDrag = fuselageDragArea / S
    return fuselageDrag


# interferance and corrections
def interference_Drag():
    viscousInter = InterferenceEffectAdCorrections.visoucsInterference(C_F_fuse, t_r, C_ci, Lambda_half)
    wingCorr = InterferenceEffectAdCorrections.wingBaseCorrection(C_F, c_r, D_F, C_L)
    # for straight or low bypass jet engines, 50% penality on the profile drag of the nacelle engine + pylon. No data for high bypass
    # upWash = InterferenceEffectAdCorrections.upwashCorrection(beta, sweepback_angle, aspect_ratio, d_2, C_L, C_L_alpha) trouver d_2 huh
    wingTailInter = InterferenceEffectAdCorrections.wingTailplaneInterference_vortexInduced(C_L_h, C_L, deps_dC_L,
                                                                                            aspect_ratio, S_tail)
    InterfDragArea = viscousInter + wingCorr + wingTailInter
    InterfDrag = InterfDragArea / S
    # for viscous interferences of the tailplane, we need to use groos area instead of wetted area for the profile drag of the tailplane.
    return InterfDrag


# Protuberances and surface imperfections
def imperfections_Drag():
    protrudingIntake = Protuberances_surfaceImperfections.protrudingIntake(Area)
    flaps_dragArea = 0  # flaps.png
    imperfection_dragArea = 0  # imperfections.png
    imperfectionDragArea = protrudingIntake + flaps_dragArea + imperfection_dragArea
    imperfectionDrag = imperfectionDragArea / S
    return imperfectionDrag


def total_Drag():
    wingDrag = wing_Drag()
    engineDrag = engine_Drag()
    tailDrag = tail_Drag()
    fuselageDrag = fuselage_Drag()
    InterfDrag = interference_Drag()
    imperfectionDrag = imperfections_Drag()
    totalDrag = wingDrag + engineDrag + tailDrag + fuselageDrag + InterfDrag + imperfectionDrag
    print(
        "wing drag: {}\nengine drag: {}\ntail drag: {}\nfuselage drag: {}\ninterference drag: {}\nimperfection drag: {}\ntotal drag:{}".format(
            wingDrag, engineDrag, tailDrag, fuselageDrag, InterfDrag, imperfectionDrag, totalDrag))
    return totalDrag


total_Drag()
