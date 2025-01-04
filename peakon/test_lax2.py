from sympy import *

u, ux, l, m = symbols('u ux l m')
mt, mx, uxx = symbols('mt mx uxx')
ut, utxx, uxxx = symbols('ut utxx uxxx')

A = Matrix([
	[-u*ux, -ux/l - u**2*m, -ux**2/l],
	[u, 1/l, ux/l - u**2*m],
	[u**2*l, u, u*ux]
])
Lt = Matrix([
	[0, mt, 0],
	[0, 0, mt],
	[0, 0, 0]
])
Ax = Matrix([
	[-ux**2-u*uxx										, -uxx/l - 2*u*ux*m - u**2*mx, -2*ux*uxx/l],
	[ux													, 0, uxx/l - 2*u*ux*m - u**2*mx],
	[2*u*ux*l, ux, ux**2 + u*uxx]
])

def novikov_compatability():
	L = Matrix([
		[0, m, -1/l],
		[0, 0, m],
		[-l, 0, 0]
	])
	expr = L*A - A*L + Lt - Ax
	print("m = u - uxx")
	pprint(expr)
	print("mt + mx*u**2 + 3*u*ux*m = 0")
	expr = simplify(expr.subs({m: u - uxx}))
	pprint(expr)
	print()


def HF_novikov_compatability():
	L = Matrix([
		[0, m, 0],
		[0, 0, m],
		[-l, 0, 0]
	])
	expr = L*A - A*L + Lt - Ax
	print("m = -uxx")
	pprint(expr)
	print("mt + mx*u**2 + 3*u*ux*m = 0")
	expr = simplify(expr.subs({m: -uxx}))
	pprint(expr)
	print()

print("Novikov Compatability")
novikov_compatability()

print("High frequency Novikov Compatability")
HF_novikov_compatability()