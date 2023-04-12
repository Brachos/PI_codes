import math


def BasicDrag(C_F, t, c, Lambda_half, S_h):
    C_D__S_h_basic = 2 * C_F * (1 + 2.75 * (t / c) * math.cos(Lambda_half) ** 2) * S_h  # F-63
    return C_D__S_h_basic


def DragDueToLift(C_L_h, Lambda_h, A_h, S_h):
    deltaC_D__S = .33 * C_L_h ** 2 / (math.cos(Lambda_h) ** 2 * math.pi * A_h) * S_h  # F-64
    return deltaC_D__S
