from sympy import *

u, ux, z, m = symbols('u ux z m')
mt, mx, uxx = symbols('mt mx uxx')
ut, utxx, uxxx = symbols('ut utxx uxxx')
C = symbols('C')

A = Matrix([
	[-u*ux, ux/z - u**2*m*z, ux**2],
	[u/z, -1/(z**2), -ux/z - u**2*m*z],
	[-u**2, u/z, u*ux]
])
Lt = Matrix([
	[0, z*mt, 0],
	[0, 0, z*mt],
	[0, 0, 0]
])
Ax = Matrix([
	[-ux**2-u*uxx, uxx/z - 2*u*ux*m*z - u**2*mx*z, 2*ux*uxx],
	[ux/z, 0, -uxx/z - 2*u*ux*m*z - u**2*mx*z],
	[-2*u*ux, ux/z, ux**2 + u*uxx]
])

def novikov_compatability():
	L = Matrix([
		[0, z*m, 1],
		[0, 0, z*m],
		[1, 0, 0]
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
		[0, z*m, 0],
		[0, 0, z*m],
		[1, 0, 0]
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