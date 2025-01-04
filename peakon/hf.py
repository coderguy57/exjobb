from sympy import *

N = 3

m = []
for i in range(N):
	m.append(Symbol('m' + str(i + 1)))

E = {}
for i in range(1,N+1):
	for j in range(1,N+1):
		if i != j:
			if (j, i) in E:
				E[(i, j)] = E[(j, i)]
			else:
				E[(i, j)] = Symbol('E' + str(i) + str(j))
		else:
			E[(i, j)] = 0

F = Matrix(N, N, lambda i, j: E[(i+1, j+1)])
P = diag(*m)
T = 2*ones(N).lower_triangular() - eye(N)

M = P*F*P

total = 0
for term in M.expand():
	total = total + term
pprint(total)

H2 = 0
for i in range(N):
	for j in range(N):
		H2 = H2 + M.minor(i, j)


def sub(expr):
	x1, x2, x3 = symbols('x1 x2 x3')
	return expr.subs(E[(1,2)], (x2-x1)).subs(E[(1,3)], (x3-x1)).subs(E[(2,3)], (x3-x2))

print(sub(H2))
pprint(M.det().simplify())

# m1, m2, m3 = m
# h2 = (1 - E[(1,2)]**2)*(m1*m2)**2 + (1 - E[(1,3)]**2)*(m1*m3)**2 + (1 - E[(2,3)]**2)*(m2*m3)**2 + \
#  	2*(E[(2,3)] - E[(1,2)]*E[(1,3)])*m1**2*m2*m3 + 2*(E[(1,2)] - E[(2,3)]*E[(1,3)])*m1*m2*m3**2
# x1, x2, x3 = symbols('x1 x2 x3')
# print(H2.subs(E[(1,2)], (x2-x1)).subs(E[(1,3)], (x3-x1)).subs(E[(2,3)], (x3-x2)))

pprint(P*F*P)