# must be considered as minimum effects
import math


def wingLiftCarryOver(eta_f, lambd, C_L_0, A):
    """
    for mid wing
    :param eta_f:
    :param lambd:
    :param C_L_0:
    :param A:
    :return:
    """
    C_D_V = (.55 * eta_f) / (1 + lambd) * (2 - math.pi * eta_f) * C_L_0 ** 2 / (math.pi * A)  # F-65
    return C_D_V


def visoucsInterference(C_F, t_r, C_ci):
    C_D__S = 1.5 * C_F * t_r * C_ci * math.cos(Lambda_half) ** 2  # F-66
    return C_D__S


def wingBaseCorrection(C_F, c_r, D_F):
    C_D__S = (-.81) * C_F * C_L * c_r * D_F  # F-67
    return C_D__S


def upwashCorrection(beta, Lambda, A, d_2, C_L, C_L_alpha):
    C_D__S = 2 * beta * math.cos(Lambda) / A * d_2 * C_L / C_L_alpha  # F-69
    return C_D__S


def wingTailplaneInterference_vortexInduced(C_L_h, C_L, deps_dC_L, A, S_h):
    C_D__S = C_L_h * C_L * (deps_dC_L - 2 / (math.pi * A)) * S_h  # F-71 can be negative
    return C_D__S


def viscousInterf():
    return 0
# NO NEED IF GROSS AREA USED IN SECTION F-3.6 INSTEAD OF WETTED AREA :DDDDDD
