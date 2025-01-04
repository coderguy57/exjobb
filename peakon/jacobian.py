from sympy import *

# Define symbols for coordinates and momenta
m1, m2, m3, x1, x2, x3 = symbols('m1 m2 m3 x1 x2 x3')

# Define the constants of motion
M1 = m1*m2*(x2 - x1) + m1*m3*(x3 - x1) + m2*m3*(x3 - x2)
M2 = (m1*m2*m3)**2*(x2 - x1)*(x3 - x1)*(x3 - x2)

A = m1**2*M1 + 2*m1**2*m2*m3*(x3-x2)
B = m2**2*M1 + 2*m2**2*m1*m3*(x3-x1)
C = m3**2*M1 + 2*m3**2*m1*m2*(x2-x1)
M3 = A*x1**0 + B*x2**0 + C*x3**0
M4 = A*x1**1 + B*x2**1 + C*x3**1
M5 = A*x1**2 + B*x2**2 + C*x3**2

M = m1*x1 + m2*x2 + m3*x3
N = m1 + m2 + m3

# Construct the functions as a list
constants_of_motion = [M1]

# Construct the Jacobian matrix
variables = [m1, m2, m3, x1, x2, x3]
Jacobian = Matrix(constants_of_motion).jacobian(variables)

print("Calculating the rank of the Jacobian matrix")
# Calculate the rank of the Jacobian
rank = Jacobian.rank()

# The Jacobian matrix and its rank
print("Jacobian matrix:")
pprint(Jacobian)
print("Rank of the Jacobian matrix: ", rank)

# pprint((M3 * M5 - M4**2 - M1**4).expand().simplify())
# pprint((M3 * M5 - M4**2).expa
# pprint((M1*M1*M1*M1).expand())
# pprint((M3 * M5 - M4**2 - M1**4).expand().simplify())

test = (M3 * M5 - M4**2 - M1**4)/((m1*m2*m3)**2)
test = test.expand().simplify()
a, b, c = symbols('a b c')
test = test.subs({m1*m2: a, m2*m3: b, m3*m1: c})
test = test.collect([a, b, c])
pprint(test)
pprint((M1*M1).expand().simplify())