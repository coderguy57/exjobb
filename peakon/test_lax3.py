from sympy import *

u, ux, l, m = symbols('u ux l m')
mt, mx, uxx = symbols('mt mx uxx')
ut, utxx, uxxx = symbols('ut utxx uxxx')

A = Matrix([
	[-u*ux, -ux/l, -ux**2/l],
	[-u*uxx*m, 1/l, ux/l],
	[u**2*l, u, u*ux]
])
Lt = Matrix([
	[0, 0, 0],
	[0, 0, 0],
	[0, 0, 0]
])
Ax = Matrix([
	[-ux**2, 0, 0],
	[ux, 0, 0],
	[2*u*ux*l, ux, ux**2]
])

def novikov_compatability():
	L = Matrix([
		[0, 0, -1/l],
		[0, 0, 0],
		[-l, 0, 0]
	])
	expr = L*A - A*L + Lt - Ax
	print("m = u - uxx")
	pprint(expr)
	print("mt - mx*u**2 + 3*u*ux*m = 0")
	expr = simplify(expr.subs({m: u - uxx}))
	pprint(expr)
	print()


def HF_novikov_compatability():
	L = Matrix([
		[0, 0, 0],
		[0, 0, 0],
		[-l, 0, 0]
	])
	expr = L*A - A*L + Lt - Ax
	print("m = -uxx")
	pprint(expr)
	print("mt - mx*u**2 + 3*u*ux*m = 0")
	expr = simplify(expr.subs({m: -uxx}))
	pprint(expr)
	print()

print("Novikov Compatability")
novikov_compatability()

print("High frequency Novikov Compatability")
HF_novikov_compatability()