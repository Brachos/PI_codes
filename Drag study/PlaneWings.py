# beta = sqrt(1-M^2): Prandtl-Glauert compressibility factor
# A: aspect ratio
# Lambda_(1/4): sweepback angle of quarter-chord line
# lambda: taper ratio
# Lambda_beta = atg(tg(Lambda)/beta): corrected sweepback angle
import math


# -------------------Method A-------------------#
# Accurate if (beta * A) < 10


def untwistedPlaneWings(C_1, C_2, C_3, taper_ratio, swept_angle, C_L, A):
    """
    vortex induced drag from wing F-11
    :param C_: Values of C_1,2,3 can be found in the fig "C_values.png" in this folder
    :param taper_ratio:
    :param swept_angle:
    :param C_L:
    :param A: aspect ratio
    :return: drag due to the wing, without the twist
    """

    eta_cp = C_1 * (1 + 2 * taper_ratio) / (3 * (1 + taper_ratio)) + (C_2 + C_3) * 4 / (
            3 * math.pi) + .001 * swept_angle * C_3
    delta = 46.264 * (eta_cp - 4 / (3 * math.pi)) ** 2
    C_D_v = (1 + delta) * C_L ** 2 / (math.pi * A)
    A = 0
    B = 0
    C = C_D_v / C_L ** 2
    return C_D_v


def wingSection(c_l_max, t, c, c_l, c_l_i, C_F, shape, Lambda_half):
    ref = 0.1 * c_l_max - .0046 * (1 + 2.75 * t / c)
    Deltac_d = 0.75 * ref * ((c_l - c_l_i) / (c_l_max - c_l_i)) ** 2
    c_d_p_min = 2 * C_F * (1 + shape * math.cos(Lambda_half) ** 2)
    c_d_p = Deltac_d + c_d_p_min
    return c_d_p


def wing3D(C_L_max, t, c, C_L, C_L_i, C_F, shape, Lambda_half, S, S_net):
    """
    This is the second simplier method

    :param C_L_max:
    :param t:
    :param c:
    :param C_L:
    :param C_L_i:
    :param C_F:
    :param shape:
    :param Lambda_half:
    :param S:
    :param S_net:
    :return:
    """
    C_D_p = wingSection(C_L_max, t, c, C_L, C_L_i, C_F * (S_net / S), shape, Lambda_half)
    return C_D_p


def wingTwist(epsilon_t, C_o_l, C_L, C_l_l):
    """
    beta * A need to be computed to find the factors
    :param epsilon_t: angle of twist
    :param C_o_l: factor found in figure "wingTwistFactors.png"
    :param C_L: lift coefficient
    :param C_l_l: actor found in figure "wingTwistFactors.png"
    :return: additionnal drag due to the twist
    """
    C_D_v = epsilon_t ** 2 * C_o_l + epsilon_t * C_L * C_l_l  # F-21
    A = 0
    B = C_D_v / C_L
    C = 0
    return C_D_v
