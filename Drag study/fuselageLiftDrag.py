import math

import computing


def VortexDragByFuselageLift(C_L, C_L_0, C_L_alpha, V_f):
    """
    If the fuselage were squared, the drag would be double. Here, with our lifting body, the drag should be less.
    :param C_L: lift ceoff
    :param C_L_0: C_L at AOA of 0 (relative to fuselage)
    :param C_L_alpha: lift-curve slope of aircraft (rad^-1)
    :param V_f: volume of the fuselage
    :return:
    """
    alpha_f = (C_L - C_L_0) / C_L_alpha
    C_DS = .15 * alpha_f ** 2 * V_f ** (2 / 3)
    return C_DS


def fuselag(C_F, S_f_wet, A_c, l_n, l_c, C_D_beta):
    phi_f = computing.phi_f(A_c, l_n, l_c)
    C_D__S_F = C_F * S_f_wet * (1 + phi_f / 2)
    D_f_eff = computing.D_f_eff(A_c)
    C_D__S = C_D__S_F + C_D_beta * math.pi / 4 * D_f_eff ** 2
    return C_D__S
