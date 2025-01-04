from sympy import symbols, Eq, solve, simplify

# Define the symbols
x, y, z, A, B, C = symbols('x y z A B C')

# Coefficients in the equations (replace these with your actual coefficients)
m, l, k = symbols('m l k')

# Define the equations
eq1 = Eq(x - A, m*(y + B))
eq2 = Eq(y - B, m*(-l*k*(x + A) + (z + C)))
eq3 = Eq(z - C, l*k*(x - A))

# Solve the system of equations
solution = solve((eq1, eq2, eq3), (x, y, z))
solution = simplify(solution)

print(solution)
