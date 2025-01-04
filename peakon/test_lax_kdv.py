from sympy import *

u, ux, z, m = symbols('u ux z m')
mt, mx, uxx = symbols('mt mx uxx')
ut, utxx, uxxx = symbols('ut utxx uxxx')
C = symbols('C')

A = Matrix([
	[-ux, 2*u+4*z],
	[2*u**2-uxx+2*u*z-4*z**2, ux]
])
Lt = Matrix([
	[0, 0],
	[ut, 0]
])
Ax = Matrix([
	[-uxx, 2*ux],
	[4*u*ux - uxxx + 2*ux*z, uxx]
])

def novikov_compatability():
	L = Matrix([
		[0, 1],
		[u-z, 0]
	])
	expr = L*A - A*L + Lt - Ax
	print("m = u - uxx")
	pprint(simplify(expr))
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
	pprint(simplify(expand(expr)))
	print("mt - mx*u**2 + 3*u*ux*m = 0")
	expr = simplify(expr.subs({m: -uxx}))
	pprint(expr)
	print()

print("Novikov Compatability")
novikov_compatability()

print("High frequency Novikov Compatability")
HF_novikov_compatability()