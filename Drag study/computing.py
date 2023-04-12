import math


def D_f_eff(A_c):
    return math.sqrt(4 / math.pi * A_c)


def Re(V, l, rho, mu):
    return V * l * mu / rho


def phi_f(A_c, l_n, l_c):
    D_f_eff = D_f_eff(A_c)
    lambda_eff = min((2 * l_n + l_c) / D_f_eff, 2 * (1 + l_n / D_f_eff))
    phi = 2.2 / lambda_eff ** (1.5) + 3.8 / lambda_eff ** 3
    return phi


def C_values_abciss(A, c_l_alpha, Lambda_quarter):
    return math.pi * 2 * A / (c_l_alpha * math.cos(Lambda_quarter))


def plug_wetted_area(l_p, D_p):
    return 0.7 * math.pi * l_p * D_p


A =  # aspect ratio A= b^2/S
c_l_alpha =  # lift curve slope of an airfoil section
Lambda_quarter =  # sweep angle of quarter-chord line, defined in A-3.1
C_abciss = C_values_abciss(A, c_l_alpha, Lambda_quarter)
print("C parameters abciss : {}".format(C_abciss))

rho =  # property of air
mu =  # property of air
v_engine =  # speed associated to engine for Re
l_engine =  #
x_transition_point_engine =  # transition point on the enigne

C_F_transitionLoc_engine = x_transition_point_engine / l_engine
C_F_Re_engine = Re(v_engine, l_engine, rho, mu)

print("Reynolds engine : {}, transition eninge: {}".format(C_F_Re_engine, C_F_transitionLoc_engine))

v_wing =  # speed associated to wing for Re
l_wing =  #
x_transition_point_wing =  # transition point on the wing

C_F_transitionLoc_wing = x_transition_point_wing / l_wing
C_F_Re_wing = Re(v_wing, l_wing, rho, mu)

print("Reynolds wing : {}, transition wing: {}".format(C_F_Re_wing, C_F_transitionLoc_wing))

v_fuselage =  # speed associated to fuselage for Re
l_fuselage =  #
x_transition_point_fuselage =  # transition point on the fuselage

C_F_transitionLoc_fuselage = x_transition_point_fuselage / l_fuselage
C_F_Re_fuselage = Re(v_fuselage, l_fuselage, rho, mu)

print("Reynolds fuselage : {}, transition fuselage: {}".format(C_F_Re_fuselage, C_F_transitionLoc_fuselage))


v_tail =  # speed associated to tail for Re
l_tail =  #
x_transition_point_tail =  # transition point on the tail

C_F_transitiontail = x_transition_point_tail / l_tail
C_F_Re_tail = Re(v_tail, l_tail, rho, mu)

print("Reynolds tail : {}, transition tail: {}".format(C_F_Re_tail, C_F_transitionLoc_tail))