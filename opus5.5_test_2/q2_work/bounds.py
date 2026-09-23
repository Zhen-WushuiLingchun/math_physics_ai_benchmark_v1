"""Rigorous finite-N bounds on E F_rec from the exact second moments."""
import numpy as np
from fractions import Fraction as Fr
from partB_moments import formulas
from common import *

def decoupling_bound(m, n, k):
    """E delta <= 1/2 sum_a sqrt(k n_a E Tr C_a^2);  E F_rec >= (1 - E delta)^2 (if E delta<=1)"""
    f = formulas(m, n, k)
    Ed = 0.5 * sum(np.sqrt(k * n[a] * float(f[a]['TrC2'])) for a in range(2))
    return Ed, max(0.0, 1 - Ed) ** 2

def collision_bound(m, n, k):
    """E F_rec >= sum_a (E q_a)^3 / (k n_a E Tr (rho_QB^(a))^2)"""
    f = formulas(m, n, k)
    return sum(float(f[a]['q']) ** 3 / (k * n[a] * float(f[a]['TrrhoQB2'])) for a in range(2))
