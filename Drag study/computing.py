import math


def D_f_eff(A_c):
    return math.sqrt(4 / math.pi * A_c)


def Re(V, l, rho, mu):
    return V * l * rho / mu


def phi_f(A_c, l_n, l_c):
    D_f_eff_ = D_f_eff(A_c)
    lambda_eff = min((2 * l_n + l_c) / D_f_eff_, 2 * (1 + l_n / D_f_eff_))
    phi = 2.2 / lambda_eff ** (1.5) + 3.8 / lambda_eff ** 3
    return phi


def C_values_abciss(A, c_l_alpha, Lambda_quarter):
    return math.pi * 2 * A / (c_l_alpha * math.cos(Lambda_quarter))


def plug_wetted_area(l_p, D_p):
    return 0.7 * math.pi * l_p * D_p


"""--------helping functions--------"""
rho = 0.364  # property of air
mu = 4.098E-5  # property of air


def C_params():
    A = 7  # aspect ratio A= b^2/S
    c_l_alpha = 0.085 / 180 * math.pi  # lift curve slope of an airfoil section
    Lambda_quarter = 36.3 / 180 * math.pi  # sweep angle of quarter-chord line, defined in A-3.1
    C_abciss = C_values_abciss(A, c_l_alpha, Lambda_quarter)
    print("C parameters abciss : {}".format(C_abciss))


def shape_engine():
    v_engine = 236.1  # speed associated to engine for Re
    l_engine = 1  #
    x_transition_point_engine = 0.01  # transition point on the enigne

    C_F_transitionLoc_engine = x_transition_point_engine / l_engine
    C_F_Re_engine = Re(v_engine, l_engine, rho, mu)

    print("Reynolds engine : {}, transition eninge: {}".format(C_F_Re_engine, C_F_transitionLoc_engine))


def shape_wing():
    v_wing = 236.1  # speed associated to wing for Re
    l_wing = 1.126  #
    x_transition_point_wing = 0.45  # transition point on the wing

    C_F_transitionLoc_wing = x_transition_point_wing / l_wing
    C_F_Re_wing = Re(v_wing, l_wing, rho, mu)

    print("Reynolds wing : {}, transition wing: {}".format(C_F_Re_wing, C_F_transitionLoc_wing))


def shape_fuselage():
    v_fuselage = 236.1  # speed associated to fuselage for Re
    l_fuselage = 5.45  #

    C_F_transitionLoc_fuselage = 0.2
    C_F_Re_fuselage = Re(v_fuselage, l_fuselage, rho, mu)

    print("Reynolds fuselage : {}, transition fuselage: {}".format(C_F_Re_fuselage, C_F_transitionLoc_fuselage))


def shape_tail():
    v_tail = 236.1  # speed associated to tail for Re
    l_tail = 0.8  #
    x_transition_point_tail = 1  # transition point on the tail

    C_F_transitionLoc_tail = .45
    C_F_Re_tail = Re(v_tail, l_tail, rho, mu)

    print("Reynolds tail : {}, transition tail: {}".format(C_F_Re_tail, C_F_transitionLoc_tail))
