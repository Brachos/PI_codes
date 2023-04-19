import math


def BasicDrag(C_F, t, c, Lambda_half, S_h):
    C_D__S_h_basic = 2 * C_F * (1 + 2.75 * (t / c) * math.cos(
        Lambda_half) ** 2) * S_h  # F-63   The reference length for Rand tic is the mean geometric tailplane chord, and it isusual to assume transition at the leadingedge.
    return C_D__S_h_basic


def DragDueToLift(C_L_h, Lambda_h, A_h, S_h):
    deltaC_D__S = .33 * C_L_h ** 2 / (math.cos(Lambda_h) ** 2 * math.pi * A_h) * S_h  # F-64, F-2.6 for C_L_h
    return deltaC_D__S


def tailVortexDrag(C_L_h, S_h, A_h):
    C_D__S = 1.02 * C_L_h ** 2 * S_h / (math.pi * A_h)
    return C_D__S
