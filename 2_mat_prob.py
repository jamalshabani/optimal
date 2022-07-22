from firedrake import *
from petsc4py import PETSc
import time
import numpy as np

start = time.time()

# Import "gmesh" mesh
mesh = Mesh("main.msh")
Id = Identity(mesh.geometric_dimension()) # Identity tensor

# Define the function spaces
V = FunctionSpace(mesh, 'CG', 1)
VV = VectorFunctionSpace(mesh, 'CG', 1)

# Create initial design
###### Begin Initial Design #####
mesh_coordinates = mesh.coordinates.dat.data[:]

# Define the points
r1 = [[0, 0], [0.2, 0], [0.4, 0], [0.6, 0], [0.8, 0], [1.0, 0]]

r2 = [[0.1, 1/12], [0.3, 1/12], [0.5, 1/12], [0.7, 1/12], [0.9, 1/12]]

r3 = [[0, 2/12], [0.2, 2/12], [0.4, 2/12], [0.6, 2/12], [0.8, 2/12], [1.0, 2/12]]

r4 = [[0.1, 3/12], [0.3, 3/12], [0.5, 3/12], [0.7, 3/12], [0.9, 3/12]]

r5 = [[0, 4/12], [0.2, 4/12], [0.4, 4/12], [0.6, 4/12], [0.8, 4/12], [1.0, 4/12]]

r = r1 + r2 + r3 + r4 + r5

def dist(x, y, a, b):
	return sqrt((x - a)**2 + (y - b)**2)

def g(s):
	if (s < 0.02):
		return 1
	if (0.02 < s < 0.03):
		return (cos((s - 0.02)*pi/0.01) + 1)/2
	if (s > 0.03):
		return 0

M = len(mesh_coordinates)

rho_array = np.zeros(M)

for i in range(len(r)):
	for j in range(M):
		x = mesh_coordinates[j][0]
		y = mesh_coordinates[j][1]
		temp = g(dist(x, y, r[i][0], r[i][1]))
		if (temp > 0):
			rho_array[j] = temp

rho = Function(V)
# rho = Constant(0.5)
rho.dat.data[:] = rho_array
x, y = SpatialCoordinate(mesh)
rho = interpolate(rho, V)

###### End Initial Design #####


# Define the constant parameters used in the problem
c = 1.0e-2
lagrange = 5.0 # Lagrange multiplier for Volume constraint
delta = 1.0e-3
epsilon = 5.0e-3

c_d_e = c / epsilon
c_m_e = c * epsilon

f = Constant((0, -1))

# Young's modulus of the beam and poisson ratio
E = 1.0
nu = 0.3 #nu poisson ratio

mu = E/(2 * (1 + nu))
_lambda = (E * nu)/((1 + nu) * (1 - 2 * nu))

# Define h(x)=x^2
def h(rho):
	return delta * (1 - rho)**2 + rho**2

# Define W(a) function i.e W(x) = x(1 - x)
def W(rho):
	return rho * (1 - rho)

# Define stress and strain tensors
def epsilon(u):
    return 0.5 * (grad(u) + grad(u).T)

def sigma(u, Id):
    return _lambda * tr(epsilon(u)) * Id + 2 * mu * epsilon(u)

# Define test function and beam displacement
v = TestFunction(VV)
u = Function(VV)

# The left side of the beam is clamped
bcs = DirichletBC(VV, Constant((0, 0)), 7)

# Define the objective function
func1 = inner(f, u) * ds(8)
func2 = c_d_e * W(rho) * dx
func3 = c_m_e * inner(grad(rho), grad(rho)) * dx
func4 = lagrange * rho * dx  # Volume constraint

J = func1 + func2 + func3 + func4

# Define the weak form for forward PDE
a_forward = h(rho) * inner(sigma(u, Id), epsilon(v)) * dx
L_forward = inner(f, v) * ds(8)
R_fwd = a_forward - L_forward

# Define the Lagrangian
a_legrange = h(rho) * inner(sigma(u, Id), epsilon(u)) * dx
L_legrange = inner(f, u) * ds(8)
R_legrange = a_legrange - L_legrange
L = J - R_legrange


def FormObjectiveGradient(tao, x, G):

	i = tao.getIterationNumber()
	if (i%10) == 0:
		File("output1/rho-{}.pvd".format(i)).write(rho)

	with rho.dat.vec as rho_vec:
		rho_vec.set(0.0)
		rho_vec.axpy(1.0, x)

	volume = assemble(rho * dx) * 3
	print("The volume fraction is {}".format(volume))

	# Solve forward PDE
	solve(R_fwd == 0, u, bcs = bcs)

	dJdrho = assemble(derivative(L, rho))
	with dJdrho.dat.vec as dJdrho_vec:
		G.set(0.0)
		G.axpy(1.0, dJdrho_vec)

	f_val = assemble(J)
	return f_val

# Setting lower and upper bounds
lb = Constant(0)
ub = Constant(1)
lb = interpolate(lb, V)
ub = interpolate(ub, V)

with lb.dat.vec as lb_vec:
	rho_lb = lb_vec

with ub.dat.vec as ub_vec:
	rho_ub = ub_vec

# Setting TAO solver
tao = PETSc.TAO().create(PETSc.COMM_SELF)
#tao.setType('bncg')
tao.setType('blmvm')
tao.setObjectiveGradient(FormObjectiveGradient, None)
tao.setVariableBounds(rho_lb, rho_ub)
tao.setFromOptions()

# Initial design guess
with rho.dat.vec as rho_vec:
	x = rho_vec.copy()

# Solve the optimization problem
tao.solve(x)
tao.destroy()

# Recover the final solution
with rho.dat.vec as rho_vec:
	rho_vec = x.copy()

# Save the final solution
File("output1/final-rho.pvd").write(rho)

end = time.time()
print("\nExecution time (in seconds):", (end - start))
