# Requires: sympy >= 1.12
import sympy as sp

# ---- helpers & basic symbols ----
n = 4
t = sp.symbols('t')

# base (time-free) symbols
x = list(sp.symbols(' '.join(f'x{i}' for i in range(1, n+1))))
m = list(sp.symbols(' '.join(f'm{i}' for i in range(1, n+1))))

# time-dependent versions
X = [sp.Function(f'x{i}')(t) for i in range(1, n+1)]
M = [sp.Function(f'm{i}')(t) for i in range(1, n+1)]

# sign on integer indices
def sgn(i, j):
    return 1 if i > j else (-1 if i < j else 0)

# u[k] and ux[k], 1-based external API
def u(k):
    k -= 1  # convert to 0-based
    return sp.Add(*(
        m[j] * (x[j] - x[k]) * sgn(j+1, k+1)
        for j in range(n)
    ))

def ux(k):
    k -= 1
    return sp.Add(*(
        m[j] * sgn(k+1, j+1)
        for j in range(n)
    ))

# time/untime substitution maps (mirror Mathematica rules)
time_subs = {x[i]: X[i] for i in range(n)} | {m[i]: M[i] for i in range(n)}
untime_subs = (
    {sp.diff(X[i], t): u(i+1)**2 for i in range(n)} |
    {sp.diff(M[i], t): -m[i]*u(i+1)*ux(i+1) for i in range(n)} |
    {X[i]: x[i] for i in range(n)} |
    {M[i]: m[i] for i in range(n)}
)

# dot[f_] := D[f /. time, t] /. untime
def dot(f):
    return sp.expand(sp.diff(sp.expand(f).xreplace(time_subs), t).xreplace(untime_subs))

# ---- build candidate conserved quantity "cand" ----
# a[j1,j2,j3,i1,i2,i3,i4] as distinct symbols for ordered index tuples
a_syms = {}
def A(j1,j2,j3,i1,i2,i3,i4):
    key = (j1,j2,j3,i1,i2,i3,i4)
    if key not in a_syms:
        a_syms[key] = sp.symbols(f"a_{j1}{j2}{j3}{i1}{i2}{i3}{i4}")
    return a_syms[key]

cand_terms = []
for i1 in range(1, n+1):
    for i2 in range(i1, n+1):
        for i3 in range(i2, n+1):
            for i4 in range(i3, n+1):
                pref = m[i1-1]*m[i2-1]*m[i3-1]*m[i4-1]
                inner = []
                for j1 in range(1, n+1):
                    for j2 in range(j1, n+1):
                        for j3 in range(j2, n+1):
                            inner.append(
                                A(j1,j2,j3,i1,i2,i3,i4)*x[j1-1]*x[j2-1]*x[j3-1]
                            )
                cand_terms.append(pref * sp.Add(*inner))
cand = sp.expand(sp.Add(*cand_terms))

# ---- build linear system: coefficients of dot(cand) vanish ----
vars_for_poly = m + x  # order matches Mathematica's Join[m /@ Range[n], x /@ Range[n]]
dot_cand = sp.expand(dot(cand))
poly = sp.Poly(dot_cand, *vars_for_poly)

# Nonzero coefficients (SymPy omits zeros); set each to 0
coeffs = poly.coeffs()
eqs = [sp.Eq(c, 0) for c in coeffs]

# Unknowns are the a[...] symbols used in cand
unknowns = list(a_syms.values())

# Solve the linear system
# (May return multiple solutions; use the first one, analogous to Mathematica's First@Solve)
sol_list = sp.solve(eqs, unknowns, dict=True)
if not sol_list:
    raise RuntimeError("No solution found for the coefficients a[...]")
sol = sol_list[0]

# ---- construct the normalized invariant and count terms ----
# Set a[1,1,1,1,1,1,2] -> 1 (normalization)
a_norm = A(1,1,1,1,1,1,2)
constM43 = sp.expand(cand.subs(sol).subs({a_norm: 1}))

# number of additive terms (like Mathematica's Length after Last expression)
num_terms = len(sp.Add.make_args(constM43))

# Optional: inspect nonzero raw coefficient expressions before solving (like Select[Flatten[cs], #!=0 &])
# nonzero_coeffs = [sp.simplify(c) for c in coeffs]  # for debugging / inspection

# Show/return results
print("constM43 =", constM43)
print("Number of terms =", num_terms)