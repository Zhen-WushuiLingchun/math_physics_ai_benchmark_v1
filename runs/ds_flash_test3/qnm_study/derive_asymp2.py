"""
Symbolically derive the asymptotic expansion of the Jost solution at x->+inf
for the Regge-Wheeler potential with general l:
   V = (1-2/rho)( l(l+1)/rho^2 - 6/rho^3 ),   x = rho + 2 ln(rho/2 - 1)
   psi_R = e^{i w x} F(x),  F = 1 + sum_{k<=N} G_k(log x) x^{-k}
Usage:  python derive_asymp2.py l N   ->  saves asymp_coeffs_l{l}_N{N}.pkl
"""
import sys
import sympy as sp
import pickle

l = int(sys.argv[1]) if len(sys.argv) > 1 else 2
N = int(sys.argv[2]) if len(sys.argv) > 2 else 8
K = N + 3

Lx, t, w = sp.symbols('Lx t omega')

# rho - 2 = (1/t)(1 + sum_{k>=1} B_k t^k)
B = [sp.symbols(f'b{k}') for k in range(1, K + 1)]
Uform = sum(B[k - 1] * t**k for k in range(1, K + 1))
eq = 2 + Uform / t + 2 * Lx - 2 * sp.log(2) + 2 * sp.log(1 + Uform)
eq_ser = sp.expand(sp.series(eq, t, 0, K + 1).removeO())
known = {}
for order in range(0, K + 1):
    coeff = sp.expand(eq_ser.coeff(t, order)).subs(known)
    vars_here = sorted({s for s in coeff.free_symbols if s in B}, key=lambda s: str(s))
    if vars_here:
        assert len(vars_here) == 1
        sol = sp.solve(sp.Poly(coeff, Lx), vars_here[0])
        known[vars_here[0]] = sp.simplify(sol[0])

U = (1 / t) * (1 + sum(known[B[k - 1]] * t**k for k in range(1, K + 1)))
rho = 2 + U
Vser = sp.expand(sp.series((1 - 2 / rho) * (l * (l + 1) / rho**2 - 6 / rho**3),
                           t, 0, K + 1).removeO())

def ddx(e):
    return sp.expand(t * sp.diff(e, Lx) - t**2 * sp.diff(e, t))

Gsym = {}
for k in range(1, N + 1):
    deg = k + 1
    Gsym[k] = sum(sp.symbols(f'g{k}_{j}') * Lx**j for j in range(0, deg + 2))
F = 1 + sum(Gsym[k] * t**k for k in range(1, N + 1))
eqF = sp.expand(ddx(ddx(F)) + 2 * sp.I * w * ddx(F) - Vser * F)
eqF = sp.expand(sp.series(eqF, t, 0, N + 3).removeO())
sol = {}
for order in range(0, N + 3):
    coeff = sp.expand(eqF.coeff(t, order)).subs(sol)
    vars_here = sorted({s for s in coeff.free_symbols if str(s).startswith('g')},
                       key=lambda s: str(s))
    if vars_here:
        eqs = sp.Poly(coeff, Lx).all_coeffs()
        s = sp.solve(eqs, vars_here, dict=True)
        assert len(s) == 1, (order, vars_here)
        for v_, val in s[0].items():
            sol[v_] = sp.cancel(sp.together(val))

Gsol = {k: sp.expand(Gsym[k].subs(sol)) for k in range(1, N + 1)}
fname = f'asymp_coeffs_l{l}_N{N}.pkl'
with open(fname, 'wb') as f:
    pickle.dump({'G': Gsol, 'N': N, 'l': l}, f)
print("saved", fname)
print("G1 =", Gsol[1])
print("G2 =", sp.simplify(Gsol[2]))
print("G3 =", sp.simplify(Gsol[3]))
