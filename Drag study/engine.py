import math


def FanCowlingDragArea_detailed(C_F, beta, S_n_wet, C_D_beta):
    factor = (1 - 31.55 * w_a_dot * theta ** (1 / 2) / (p_inf * M_inf * D_h ** 2)) / (D_n - D_h) ** 2 - 1
    phi_n = 0.33 * (D_n - D_h) / (beta * l_n) * (1 + 1.75 * factor)
    C_D__S = C_F * (beta * (1 + phi_n) ** 5 / 3 + (1 - beta)) * S_n_wet + C_D_beta * math.pi / 4 * D_n ** 2


def FanCowlingDragArea_simple(C_F, S_n_wet):
    C_D__S = 1.25 * C_F * S_n_wet  # F-58
    return C_D__S


def gasGeneratorCowlingDragComponent_ManInfo(C_F, M_f, M_inf, S_g_wet, C_D_beta, D_g):
    C_D__S = C_F * (M_f / M_inf) ** (11 / 6) * ((1 + .116 * M_inf ** 2) / (1 + .116 * M_f ** 2)) ** (
            2 / 3) * S_g_wet + C_D_beta * math.pi / 4 * D_g ** 2  # F-59
    return C_D__S


def plugDragArea(C_F, S_p_wet, M_g, M_inf):
    C_D__S = C_F * S_p_wet * (M_g / M_inf) ** (11 / 6) * ((1 + .116 * M_inf ** 2) / (1 + .116 * M_g ** 2)) ** (
            2 / 3)  # F-60
    return C_D__S


def pylonDrag(C_F, t, c, Lambda_py, S_py_wet):
    C_D__S = C_F * (1 + 2.75 * t / c * math.cos(Lambda_py) ** 2) * S_py_wet  # F-61
    return C_D__S
