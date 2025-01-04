from sympy import *

u, ux, l, m, x = symbols('u ux l m x')
mt, mx, uxx = symbols('mt mx uxx')
ut, utxx, uxxx = symbols('ut utxx uxxx')

A = Matrix([
	[ux*(ux*x-u), -ux/(l*x)-u**2*m/x, -ux**2/l],
	[-(ux*x-u)*x + l*u**2*m*x**2, 1/l, ux*x/l-u**2*m*x],
	[l*(ux*x-u)*(ux*x-u), -ux+u/x-l*u**2*m, -ux*(ux*x-u)]
])
Lt = Matrix([
	[0, mt/x, 0],
	[-l*mt*x**2, 0, mt*x],
	[0, l*mt, 0]
])
Ax = Matrix([
	[uxx*(ux*x-u)+ux*uxx*x, -uxx/(l*x)+ux/(x**2*l) - 2*u*ux*m/x - u**2*mx/x + u**2*m/x**2, -2*ux*uxx/l],
	[-uxx*x**2-(ux*x-u) + l*(2*u*ux*m*x**2 + u**2*mx*x**2 + 2*u**2*m*x), 0, uxx*x/l + ux/l - 2*u*ux*m*x - u**2*mx*x - u**2*m],
	[2*l*(ux*x-u)*(uxx*x), -uxx +ux/x-u/x**2 -l*(2*u*ux*m + u**2*mx), -uxx*(ux*x-u)-ux*uxx*x]
])

def HF_novikov_compatability():
	L = Matrix([
		[0, m/x, 0],
		[-l*m*x**2, 1/x, m*x],
		[0, l*m, 0]
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