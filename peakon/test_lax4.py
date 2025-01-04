from sympy import *

u, ux, l, m, x = symbols('u ux l m x')
mt, mx, uxx = symbols('mt mx uxx')
ut, utxx, uxxx = symbols('ut utxx uxxx')

A = Matrix([
	[ux*(ux*x-u), -ux/l-u**2*m, -ux**2/l],
	[-(ux*x-u) + l*u**2*m*x, 1/l, ux/l-u**2*m],
	[l*(ux*x-u)*(ux*x-u), -(ux*x-u)-l*u**2*m*x, -ux*(ux*x-u)]
])
Lt = Matrix([
	[0, mt, 0],
	[-l*mt*x, 0, mt],
	[0, l*mt*x, 0]
])
Ax = Matrix([
	[uxx*(ux*x-u)+ux*uxx*x, -uxx/l - 2*u*ux*m - u**2*mx, -2*ux*uxx/l],
	[-uxx*x + l*(2*u*ux*m*x + u**2*mx*x + u**2*m), 0, uxx/l - 2*u*ux*m - u**2*mx],
	[2*l*(ux*x-u)*(uxx*x), -uxx*x-l*(2*u*ux*m*x + u**2*mx*x + u**2*m), -uxx*(ux*x-u)-ux*uxx*x]
])

def HF_novikov_compatability():
	L = Matrix([
		[0, m, 0],
		[-l*m*x, 0, m],
		[0, l*m*x, 0]
	])
	expr = L*A - A*L + Lt - Ax
	print("m = -uxx")
	pprint(expr)
	print("mt + mx*u**2 + 3*u*ux*m = 0")
	expr = simplify(expr.subs({m: -uxx}))
	pprint(expr)
	print()

print("High frequency Novikov Compatability")
HF_novikov_compatability()