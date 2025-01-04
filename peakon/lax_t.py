from sympy import *

u, ux, l, m = symbols('u ux l m')
mt, mx, uxx = symbols('mt mx uxx')
ut, utxx, uxxx = symbols('ut utxx uxxx')
C = symbols('C')

A = Matrix([
	[-u*ux, -ux/l, -ux**2/l],
	[u, 1/l, ux/l],
	[u**2*l, u, u*ux]
])
p = Matrix([[0], [-ux], [1]])

pprint(A*A*A*p)
P, D = A.diagonalize()
pprint(P*(D*D)*P.inv())
pprint(A*A)