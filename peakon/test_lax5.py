from sympy import *

u, ux, l, m, x = symbols('u ux l m x')
mt, mx, uxx = symbols('mt mx uxx')
ut, utxx, uxxx = symbols('ut utxx uxxx')
s, sx = symbols('s sx')

A = Matrix([
	[ux*(ux*x-u), -ux*s/l-u**2*m*s, -ux**2/l],
	[-(ux*x-u)*s + l*u**2*m*x*s, 1/l, ux*s/l-u**2*m*s],
	[l*(ux*x-u)*(ux*x-u), -(ux*x-u)*s-l*u**2*m*x*s, -ux*(ux*x-u)]
])
Lt = Matrix([
	[0, mt*s, 0],
	[-l*mt*x*s, 0, mt*s],
	[0, l*mt*x*s, 0]
])
Ax = Matrix([
	[uxx*(ux*x-u)+ux*uxx*x, -uxx*s/l-ux*sx/l - 2*u*ux*m*s - u**2*mx*s - u**2*m*sx, -2*ux*uxx/l],
	[-uxx*x*s-(ux*x-u)*sx + l*(2*u*ux*m*x*s + u**2*mx*x*s + u**2*m*s + u**2*m*x*sx), 0, uxx*s/l + ux*sx/l - 2*u*ux*m*s - u**2*mx*s - u**2*m*sx],
	[2*l*(ux*x-u)*(uxx*x), -uxx*x*s-(ux*x-u)*sx -l*(2*u*ux*m*x*s + u**2*mx*x*s + u**2*m*s + u**2*m*x*sx), -uxx*(ux*x-u)-ux*uxx*x]
])

def HF_novikov_compatability():
	L = Matrix([
		[0, m*s, 0],
		[-l*m*x*s, s*sx, m*s],
		[0, l*m*x*s, 0]
	])
	expr = L*A - A*L + Lt - Ax
	print("m = -uxx")
	pprint(expr)
	print("mt + mx*u**2 + 3*u*ux*m = 0")
	expr = simplify(expr.subs({m: -uxx}))
	expr = simplify(expr.subs({s**2: 1}))
	# expr = simplify(expr.subs({sx: 0}))
	pprint(expr)
	print()

print("High frequency Novikov Compatability")
HF_novikov_compatability()