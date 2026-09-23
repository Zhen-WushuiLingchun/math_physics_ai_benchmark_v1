"""
Symbolically derive the asymptotic expansion of the RW Jost solution at x->+inf:
  psi_R = e^{i w x} F(x),  F(x) = 1 + sum_k G_k(log x) x^{-k}
with V(x) = 6(rho-2)(rho-1)/rho^4,  x = rho + 2 ln(rho/2 - 1).
Also produces the expansion of V itself.
"""
import sympy as sp

Lx, t, w = sp.symbols('Lx t omega')
# x = 1/t, Lx = log x.
K = 8   # order of expansion in t

# rho - 2 = (1/t) * (1 + sum_{k>=1} B_k(Lx) t^k)
# From x = rho + 2 ln(rho/2 - 1), with rho = 2 + U, U = (1/t)(1+sum B_k t^k):
#   1/t = 2 + U + 2 ln(U/2)
# Multiply by t: 1 = 2t + (1 + sum B_k t^k) + 2t ln(1/t) + 2t ln((1+sum B_k t^k)/2)
#   = 2t + 1 + sum B_k t^k - 2t Lx - 2t ln2 + 2t ln(1+sum B_k t^k)
# => 0 = 2 + sum B_k t^{k-1} - 2 Lx - 2 ln2 + 2 ln(1 + sum B_k t^k)
B = [sp.symbols(f'b{k}') for k in range(1, K+1)]
Uform = sum(B[k-1]*t**k for k in range(1, K+1))
# equation as series in t (Lx symbolic)
eq = 2 + Uform/t + 2*Lx - 2*sp.log(2) + 2*sp.log(1 + Uform)
eq_ser = sp.series(eq, t, 0, K+1).removeO()
eq_ser = sp.expand(eq_ser)
# collect order by order
known = {}
for order in range(0, K+1):
    coeff = sp.expand(eq_ser.coeff(t, order))
    coeff = coeff.subs(known)
    # solve for the unknown coefficient of t^order  (a polynomial in Lx)
    # find which b's appear with t^order
    vars_here = sorted({s for s in coeff.free_symbols if s in B}, key=lambda s: str(s))
    if vars_here:
        assert len(vars_here) == 1, (order, vars_here, coeff)
        sol = sp.solve(sp.Poly(coeff, Lx), vars_here[0])
        # coeff is linear in Lx poles? use solve on the polynomial in Lx
        # ensure identically zero as polynomial in Lx
        sol_expr = sp.simplify(sol[0])
        known[vars_here[0]] = sol_expr
    else:
        pass
print("B_k(Lx):")
for k in range(1, K+1):
    print(f"  B{k} = {sp.simplify(known[B[k-1]])}")

# now assemble rho(x) and V(x) as series in t (Lx symbolic)
U = (1/t)*(1 + sum(known[B[k-1]]*t**k for k in range(1, K+1)))
rho = 2 + U
Vser = sp.series(6*(rho-2)*(rho-1)/rho**4, t, 0, K+1).removeO()
Vser = sp.expand(Vser)
print()
print("V(x) series (coefficients of t^k):")
for k in range(2, K+1):
    print(f"  cV{k} = {sp.expand(Vser.coeff(t,k))}")

# ---------------- solution series: F = 1 + sum G_k(Lx) t^k
# d/dx = t d/dLx - t^2 d/dt
def ddx(expr):
    return sp.expand(t*sp.diff(expr, Lx) - t**2*sp.diff(expr, t))

N = 6   # solve F to order t^N
Gsym = {}
for k in range(1, N+1):
    # allow terms G_k = sum_j g[k,j] Lx^j , j = 0..k+1 (guess structure)
    deg = k + 1
    Gsym[k] = sum(sp.symbols(f'g{k}_{j}')*Lx**j for j in range(0, deg+1))
F = 1 + sum(Gsym[k]*t**k for k in range(1, N+1))
eqF = sp.expand(ddx(ddx(F)) + 2*sp.I*w*ddx(F) - Vser*F)
eqF = sp.series(eqF, t, 0, N+2).removeO()
eqF = sp.expand(eqF)
sol = {}
for order in range(0, N+2):
    coeff = sp.expand(eqF.coeff(t, order))
    coeff = coeff.subs(sol)
    vars_here = sorted({s for s in coeff.free_symbols if str(s).startswith('g')},
                       key=lambda s: str(s))
    if vars_here:
        eqs = sp.Poly(coeff, Lx).all_coeffs()
        s = sp.solve(eqs, vars_here, dict=True)
        assert len(s) == 1
        for v_, val in s[0].items():
            sol[v_] = sp.simplify(val)
print()
print("solution coefficients G_k(Lx):")
for k in range(1, N+1):
    Gk = Gsym[k].subs(sol)
    print(f"  G{k} = {sp.expand(Gk)}")

# save for later numeric use
import pickle
with open('asymp_coeffs.pkl','wb') as f:
    pickle.dump({'B': {k: known[B[k-1]] for k in range(1,K+1)},
                 'G': {k: Gsym[k].subs(sol) for k in range(1,N+1)},
                 'N': N, 'K': K}, f)
print("\nsaved to asymp_coeffs.pkl")
