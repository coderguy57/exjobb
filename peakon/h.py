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
			E[(i, j)] = 1

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

pprint(H2)
pprint(M.det().simplify())

m1, m2, m3 = m
h2 = (1 - E[(1,2)]**2)*(m1*m2)**2 + (1 - E[(1,3)]**2)*(m1*m3)**2 + (1 - E[(2,3)]**2)*(m2*m3)**2 + \
 	2*(E[(2,3)] - E[(1,2)]*E[(1,3)])*m1**2*m2*m3 + 2*(E[(1,2)] - E[(2,3)]*E[(1,3)])*m1*m2*m3**2

# pprint((H2 - h2).simplify())

pprint(P*F*P)