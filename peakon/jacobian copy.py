from sympy import symbols, Matrix, pprint

# Define symbols for coordinates and momenta
m1, m2, x1, x2 = symbols('m1 m2 x1 x2')

# Define the constants of motion
M1 = m1*m2*(x2 - x1)

A = m1**2
B = m2**2
M2 = A*x1**0 + B*x2**0
M3 = A*x1**1 + B*x2**1
M4 = A*x1**2 + B*x2**2

# Construct the functions as a list
constants_of_motion = [M1, M2, M3, M4]

# Construct the Jacobian matrix
variables = [m1, m2, x1, x2]
Jacobian = Matrix(constants_of_motion).jacobian(variables)

print("Calculating the rank of the Jacobian matrix")
# Calculate the rank of the Jacobian
rank = Jacobian.rank()

# The Jacobian matrix and its rank
print("Jacobian matrix:")
pprint(Jacobian)
print("Rank of the Jacobian matrix: ", rank)

pprint((M4*M2 - M3**2 - M1*M1).simplify())