from sympy import *

# Define symbols
m1, m2, m3, m4, x1, x2, x3, x4 = symbols('m1 m2 m3 m4 x1 x2 x3 x4')
t = symbols('t')

def sign(k1, k2):
	return 0 if k1 == k2 else (1 if k1 > k2 else -1)
def u(k):
	x = ([x1, x2, x3, x4])[k-1]
	return (m1*(x - x1)*sign(k,1) + m2*(x - x2)*sign(k,2) + m3*(x - x3)*sign(k,3) + m4*(x - x4)*sign(k,4)).simplify()
def ux(k):
	return m1*sign(k, 1) + m2*sign(k, 2) + m3*sign(k, 3) + m4*sign(k, 4)

# Define expressions for x1, x2, x3, m1, m2, m3 over time
x1_dot = u(1)**2
x2_dot = u(2)**2
x3_dot = u(3)**2
x4_dot = u(4)**2

m1_dot = -m1 * ux(1) * u(1)
m2_dot = -m2 * ux(2) * u(2)
m3_dot = -m3 * ux(3) * u(3)
m4_dot = -m4 * ux(4) * u(4)

# Define the expression whose derivative we want to find
exprs = [
-1*m3**2*m4**2*x4**2-2*m1*m3*m4**2*x4**2-1*m1**2*m4**2*x4**2+2*m3**2*m4**2*x3*x4-2*m1*m2*m4**2*x4**2+2*m2*m3*m4**2*x3*x4-2*m1*m3**2*m4*x3*x4-2*m1*m2**2*m3*x2*x3+2*m1*m2*m3**2*x1*x3+2*m1**2*m3**2*x1*x3+2*m1**2*m2*m3*x1*x3+2*m1*m2**2*m4*x2**2-2*m1**2*m2*m3*x1**2+2*m1*m2**2*m3*x2**2+2*m1**2*m3*m4*x1*x3+2*m1*m3*m4**2*x1*x4-2*m1*m2*m4**2*x1*x2+2*m2*m3**2*m4*x2*x4-1*m1**2*m2**2*x1**2-2*m2*m3*m4**2*x2*x3+2*m1*m2**2*m3*x1*x3+2*m2**2*m4**2*x2*x4+2*m1*m3*m4**2*x3*x4-2*m1*m3**2*m4*x1*x3-2*m1*m2**2*m3*x1*x2-4*m1*m2*m3*m4*x1*x2-1*m1**2*m4**2*x1**2-2*m2*m3**2*m4*x3*x4+2*m1**2*m2*m4*x1*x2-1*m2**2*m4**2*x4**2-2*m1*m2*m3**2*x3**2+2*m1**2*m2*m3*x1*x2-2*m1*m2**2*m4*x1*x2-2*m1**2*m2*m3*x2*x3+2*m1**2*m4**2*x1*x4-2*m1*m2*m3**2*x1*x2-1*m1**2*m3**2*x1**2+2*m1*m2*m3**2*x2*x3-1*m2**2*m4**2*x2**2-2*m1**2*m3*m4*x1**2-2*m1*m2**2*m4*x2*x4+2*m2**2*m3**2*x2*x3+4*m1*m2*m3*m4*x2*x3+2*m2**2*m3*m4*x2*x3-2*m2*m3**2*m4*x2*x3-2*m2*m3*m4**2*x4**2-1*m2**2*m3**2*x3**2-1*m2**2*m3**2*x2**2+2*m1*m3**2*m4*x3**2+2*m2*m3**2*m4*x3**2-1*m3**2*m4**2*x3**2-2*m1**2*m2*m4*x1**2-1*m1**2*m3**2*x3**2+2*m1*m2*m4**2*x2*x4+2*m1**2*m2*m4*x1*x4+2*m1*m2**2*m4*x1*x4+2*m1**2*m3*m4*x1*x4-2*m1*m3*m4**2*x1*x3+4*m1*m2*m3*m4*x1*x4-1*m1**2*m2**2*x2**2-2*m2**2*m3*m4*x2**2-4*m1*m2*m3*m4*x3*x4+2*m1*m3**2*m4*x1*x4+2*m1**2*m2**2*x1*x2+2*m1*m2*m4**2*x1*x4-2*m1**2*m2*m4*x2*x4+2*m2**2*m3*m4*x2*x4+2*m2*m3*m4**2*x2*x4-2*m1**2*m3*m4*x3*x4-2*m2**2*m3*m4*x3*x4,
-1*m1**2*m3**2*m4**2*x3*x4**2+1*m2**2*m3**2*m4**2*x2*x4**2+2*m1*m2*m3**2*m4**2*x2*x4**2-2*m1**2*m2*m3*m4**2*x2*x4**2-1*m1**2*m2**2*m4**2*x2*x4**2+1*m1**2*m3**2*m4**2*x1*x4**2+2*m1**2*m2*m3*m4**2*x1*x4**2+2*m1*m2*m3**2*m4**2*x3**2*x4-2*m1*m2*m3**2*m4**2*x2*x3*x4+2*m1**2*m2*m3*m4**2*x2*x3*x4-2*m1**2*m2**2*m3*m4*x2*x3*x4+2*m1*m2*m3**2*m4**2*x1*x3*x4+1*m1**2*m3**2*m4**2*x3**2*x4-2*m1**2*m2*m3*m4**2*x1*x3*x4+2*m1**2*m2**2*m3*m4*x1*x3*x4-1*m2**2*m3**2*m4**2*x2**2*x4+2*m1*m2**2*m3**2*m4*x2**2*x4-2*m1*m2*m3**2*m4**2*x1*x2*x4+2*m1*m2**2*m3**2*m4*x1*x2*x3+2*m1*m2**2*m3**2*m4*x1*x3*x4-2*m1**2*m2*m3*m4**2*x1**2*x4+2*m1**2*m2**2*m3*m4*x1*x2*x3+1*m1**2*m2**2*m3**2*x2**2*x3-2*m1*m2*m3**2*m4**2*x3*x4**2-1*m1**2*m2**2*m3**2*x1**2*x3-1*m1**2*m2**2*m4**2*x1*x2**2+1*m2**2*m3**2*m4**2*x3**2*x4-2*m1*m2**2*m3**2*m4*x2*x3*x4-2*m1**2*m2**2*m3*m4*x1*x2**2-1*m1**2*m2**2*m3**2*x1*x2**2+1*m1**2*m3**2*m4**2*x1**2*x3+2*m1**2*m2*m3*m4**2*x1**2*x3+1*m1**2*m2**2*m4**2*x1*x4**2+1*m1**2*m2**2*m4**2*x2**2*x4+2*m1**2*m2**2*m3*m4*x1**2*x2+1*m1**2*m2**2*m3**2*x1**2*x2-2*m1*m2**2*m3**2*m4*x1*x2*x4-2*m1**2*m2**2*m3*m4*x1**2*x3-2*m1*m2**2*m3**2*m4*x2**2*x3+1*m1**2*m2**2*m4**2*x1**2*x2-2*m1**2*m2*m3*m4**2*x1*x2*x3-2*m1*m2**2*m3**2*m4*x1*x3**2+2*m1*m2*m3**2*m4**2*x1*x2*x3+2*m1**2*m2**2*m3*m4*x2**2*x4+2*m1**2*m2*m3*m4**2*x1*x2*x4+1*m2**2*m3**2*m4**2*x2**2*x3+1*m1**2*m2**2*m3**2*x1*x3**2-1*m1**2*m3**2*m4**2*x1*x3**2-1*m1**2*m2**2*m3**2*x2*x3**2-1*m2**2*m3**2*m4**2*x3*x4**2-1*m1**2*m3**2*m4**2*x1**2*x4-2*m1**2*m2**2*m3*m4*x1*x2*x4-2*m1*m2*m3**2*m4**2*x1*x3**2+2*m1*m2**2*m3**2*m4*x2*x3**2-1*m2**2*m3**2*m4**2*x2*x3**2-1*m1**2*m2**2*m4**2*x1**2*x4,
-1*m1**2*m2**2*m3**2*m4**2*x2*x3*x4**2+1*m1**2*m2**2*m3**2*m4**2*x1*x3*x4**2-1*m1**2*m2**2*m3**2*m4**2*x1**2*x2*x3-1*m1**2*m2**2*m3**2*m4**2*x1*x2*x3**2+1*m1**2*m2**2*m3**2*m4**2*x1*x2**2*x3+1*m1**2*m2**2*m3**2*m4**2*x2*x3**2*x4+1*m1**2*m2**2*m3**2*m4**2*x1**2*x3**2-1*m1**2*m2**2*m3**2*m4**2*x1*x2**2*x4+2*m1**2*m2**2*m3**2*m4**2*x1*x2*x3*x4+1*m1**2*m2**2*m3**2*m4**2*x1**2*x2*x4-1*m1**2*m2**2*m3**2*m4**2*x1**2*x3*x4-1*m1**2*m2**2*m3**2*m4**2*x1*x3**2*x4-1*m1**2*m2**2*m3**2*m4**2*x2**2*x3*x4-1*m1**2*m2**2*m3**2*m4**2*x1*x2*x4**2+1*m1**2*m2**2*m3**2*m4**2*x2**2*x4**2,
-1*m3*m4*x4-1*m2*m4*x4-1*m1*m4*x4-1*m2*m3*x3+1*m3*m4*x3+1*m2*m4*x2+1*m2*m3*x2-1*m1*m3*x3-1*m1*m2*x2+1*m1*m4*x1+1*m1*m3*x1+1*m1*m2*x1
]

def test(expression):
	# Compute the total derivative w.r.t. time
	return 	diff(expression, x1)*x1_dot + diff(expression, x2)*x2_dot + diff(expression, x3)*x3_dot + diff(expression, x4)*x4_dot + \
			diff(expression, m1)*m1_dot + diff(expression, m2)*m2_dot + diff(expression, m3)*m3_dot + diff(expression, m4)*m4_dot

# for expr in exprs:
# 	print(test(expr).expand())

expr1 = exprs[0]
expr2 = exprs[1]
expr3 = exprs[2]
expr4 = exprs[3]

# pprint((expr4*expr4+expr1).expand().simplify())

# pprint((expr4*expr4).expand().simplify())
# pprint((expr2).expand().simplify())
# pprint((expr3).expand().simplify())

expr = m1*m2*(x2 - x1) * m1*m3*(x3 - x1) * m1*m4*(x4 - x1) * m2*m3*(x3 - x2) * m2*m4*(x4 - x2) * m3*m4*(x4 - x3)

# pprint(expr3.simplify())
# pprint(test(expr).simplify())

def r(a, b):
	if a > b:
		a, b = b, a
	a -= 1
	b -= 1
	m = [m1, m2, m3, m4]
	x = [x1, x2, x3, x4]
	return m[a]*m[b]*(x[b] - x[a])

# for expr in exprs:
# 	pprint(expr.expand().simplify())

t1 = r(1, 2) * r(2, 3) * r(3, 4) * r(4, 1)
t2 = r(1, 2) * r(2, 4) * r(4, 3) * r(3, 1)
t3 = r(1, 3) * r(3, 2) * r(2, 4) * r(4, 1)

# out = t1 + t2 + t3
# pprint(out.expand().simplify())
# pprint(t1.expand().simplify())
# pprint(test(t1).expand().simplify())

H1 = 2*m1*m2*(-x1 + x2) + 2*m1*m3*(-x1 + x3) + 2*m1*m4*(-x1 + x4) + 2*m2*m3*(-x2 + x3) + 2*m2*m4*(-x2 + x4) + 2*m3*m4*(-x3 + x4)
H2 = -m1**2*m2**2*(-x1 + x2)**2 - 2*m1**2*m2*m3*(-x1 + x2)*(-x1 + x3) - 2*m1**2*m2*m4*(-x1 + x2)*(-x1 + x4) - m1**2*m3**2*(-x1 + x3)**2 - 2*m1**2*m3*m4*(-x1 + x3)*(-x1 + x4) - m1**2*m4**2*(-x1 + x4)**2 + 2*m1*m2**2*m3*(-x1 + x2)*(-x2 + x3) + 2*m1*m2**2*m4*(-x1 + x2)*(-x2 + x4) - 2*m1*m2*m3**2*(-x1 + x3)*(-x2 + x3) + 4*m1*m2*m3*m4*(-x1 + x2)*(-x3 + x4) - 4*m1*m2*m3*m4*(-x1 + x4)*(-x2 + x3) - 2*m1*m2*m4**2*(-x1 + x4)*(-x2 + x4) + 2*m1*m3**2*m4*(-x1 + x3)*(-x3 + x4) - 2*m1*m3*m4**2*(-x1 + x4)*(-x3 + x4) - m2**2*m3**2*(-x2 + x3)**2 - 2*m2**2*m3*m4*(-x2 + x3)*(-x2 + x4) - m2**2*m4**2*(-x2 + x4)**2 + 2*m2*m3**2*m4*(-x2 + x3)*(-x3 + x4) - 2*m2*m3*m4**2*(-x2 + x4)*(-x3 + x4) - m3**2*m4**2*(-x3 + x4)**2
H3 = 2*m1**2*m2**2*m3**2*(-x1 + x2)*(-x1 + x3)*(-x2 + x3) - 2*m1**2*m2**2*m3*m4*(-x1 + x2)**2*(-x3 + x4) + 2*m1**2*m2**2*m3*m4*(-x1 + x2)*(-x1 + x3)*(-x2 + x4) + 2*m1**2*m2**2*m3*m4*(-x1 + x2)*(-x1 + x4)*(-x2 + x3) + 2*m1**2*m2**2*m4**2*(-x1 + x2)*(-x1 + x4)*(-x2 + x4) - 2*m1**2*m2*m3**2*m4*(-x1 + x2)*(-x1 + x3)*(-x3 + x4) + 2*m1**2*m2*m3**2*m4*(-x1 + x3)**2*(-x2 + x4) - 2*m1**2*m2*m3**2*m4*(-x1 + x3)*(-x1 + x4)*(-x2 + x3) + 2*m1**2*m2*m3*m4**2*(-x1 + x2)*(-x1 + x4)*(-x3 + x4) + 2*m1**2*m2*m3*m4**2*(-x1 + x3)*(-x1 + x4)*(-x2 + x4) - 2*m1**2*m2*m3*m4**2*(-x1 + x4)**2*(-x2 + x3) + 2*m1**2*m3**2*m4**2*(-x1 + x3)*(-x1 + x4)*(-x3 + x4) + 2*m1*m2**2*m3**2*m4*(-x1 + x2)*(-x2 + x3)*(-x3 + x4) + 2*m1*m2**2*m3**2*m4*(-x1 + x3)*(-x2 + x3)*(-x2 + x4) - 2*m1*m2**2*m3**2*m4*(-x1 + x4)*(-x2 + x3)**2 - 2*m1*m2**2*m3*m4**2*(-x1 + x2)*(-x2 + x4)*(-x3 + x4) + 2*m1*m2**2*m3*m4**2*(-x1 + x3)*(-x2 + x4)**2 - 2*m1*m2**2*m3*m4**2*(-x1 + x4)*(-x2 + x3)*(-x2 + x4) - 2*m1*m2*m3**2*m4**2*(-x1 + x2)*(-x3 + x4)**2 + 2*m1*m2*m3**2*m4**2*(-x1 + x3)*(-x2 + x4)*(-x3 + x4) + 2*m1*m2*m3**2*m4**2*(-x1 + x4)*(-x2 + x3)*(-x3 + x4) + 2*m2**2*m3**2*m4**2*(-x2 + x3)*(-x2 + x4)*(-x3 + x4)
H4 = m1**2*m2**2*m3**2*m4**2*((-x1 + x2)**2*(-x3 + x4)**2 - 2*(-x1 + x2)*(-x1 + x3)*(-x2 + x4)*(-x3 + x4) - 2*(-x1 + x2)*(-x1 + x4)*(-x2 + x3)*(-x3 + x4) + (-x1 + x3)**2*(-x2 + x4)**2 - 2*(-x1 + x3)*(-x1 + x4)*(-x2 + x3)*(-x2 + x4) + (-x1 + x4)**2*(-x2 + x3)**2)

# pprint(test(H1).expand().simplify())
# pprint(test(H2).expand().simplify())
# pprint(test(H3).expand().simplify())
# pprint(test(H4).expand().simplify())

# pprint((H1).expand().simplify())
# pprint((H2).expand().simplify())
# print((H3).expand().simplify())
# pprint((H4).expand().simplify())

pprint((2*expr2 + H3).expand().simplify())