def parse():
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('-tao_monitor', '--tao_monitor', action='store_true', help = 'TAO monitor')
	parser.add_argument('-tao_max_it', '--tao_max_it', type = int, default = 100, help = 'Number of TAO iterations')
	parser.add_argument('-l1', '--lagrange_r1', type = float, default = 1.0, help = 'Lagrange multiplier for material 1')
	parser.add_argument('-l2', '--lagrange_r2', type = float, default = 0.1, help = 'Lagrange multiplier for material 2')
	parser.add_argument('-k', '--kappa', type = float, default = 1.0e-2, help = 'Weight of Modica-Mortola')
	parser.add_argument('-e', '--epsilon', type = float, default = 5.0e-3, help = 'Phase-field regularization parameter')
	parser.add_argument('-o', '--output', type = str, default = 'output1', help = 'Output folder')
	parser.add_argument('-m', '--mesh', type = str, default = 'main.msh', help = 'Dimensions of meshed beam')
	parser.add_argument('-e1', '--er1modulus', type = float, default = 0.1, help = 'Elastic Modulus for material 1')
	parser.add_argument('-e2', '--er2modulus', type = float, default = 1.0, help = 'Elastic Modulus for material 2')
	parser.add_argument('-p', '--power_p', type = float, default = 2.0, help = 'Power for elasticity interpolation')
	options = parser.parse_args()
	return options

options = parse()

from firedrake import *
from petsc4py import PETSc
import time
import numpy as np

start = time.time()

# Import gmesh
mesh = Mesh(options.mesh)
Id = Identity(mesh.geometric_dimension()) #Identity tensor

# Define the function spaces
V = FunctionSpace(mesh, 'CG', 1)
VV = VectorFunctionSpace(mesh, 'CG', 1, dim = 2)

# Create initial design
###### Begin Initial Design #####
mesh_coordinates = mesh.coordinates.dat.data[:]

# Define the points
v0 = [[0, 0], [0.2, 0], [0.4, 0], [0.6, 0], [0.8, 0], [1.0, 0]]
r0 = [[0.1, 0], [0.3, 0], [0.5, 0], [0.7, 0], [0.9, 0]]

v1 = [[0.1, 1/12], [0.3, 1/12], [0.5, 1/12], [0.7, 1/12], [0.9, 1/12]]
r1 = [[0, 1/12], [0.2, 1/12], [0.4, 1/12], [0.6, 1/12], [0.8, 1/12], [1.0, 1/12]]

v2 = [[0, 2/12], [0.2, 2/12], [0.4, 2/12], [0.6, 2/12], [0.8, 2/12], [1.0, 2/12]]
r2 = [[0.1, 2/12], [0.3, 2/12], [0.5, 2/12], [0.7, 2/12], [0.9, 2/12]]

v3 = [[0.1, 3/12], [0.3, 3/12], [0.5, 3/12], [0.7, 3/12], [0.9, 3/12]]
r3 = [[0, 3/12], [0.2, 3/12], [0.4, 3/12], [0.6, 3/12], [0.8, 3/12], [1.0, 3/12]]

v4 = [[0, 4/12], [0.2, 4/12], [0.4, 4/12], [0.6, 4/12], [0.8, 4/12], [1.0, 4/12]]
r4 = [[0.1, 4/12], [0.3, 4/12], [0.5, 4/12], [0.7, 4/12], [0.9, 4/12]]

v = v0 + v1 + v2 + v3 + v4
r = r0 + r1 + r2 + r3 + r4
s = v + r

def dist(x, y, a, b):
	return sqrt((x - a)**2 + (y - b)**2)

def g_s(s):
	if (s < 0.02):
		return 0
	if (0.02 < s < 0.03):
		return (cos((s - 0.03)*pi/0.01) + 1)/2
	if (s > 0.03):
		return 1

def g_r(s):
	if (s < 0.02):
		return 1
	if (0.02 < s < 0.03):
		return (cos((s - 0.02)*pi/0.01) + 1)/2
	if (s > 0.03):
		return 0

M = len(mesh_coordinates)
rho2_array = np.ones(M)
rho3_array = np.zeros(M)

for i in range(len(s)):
	for j in range(M):
		x = mesh_coordinates[j][0]
		y = mesh_coordinates[j][1]
		temp = g_s(dist(x, y, s[i][0], s[i][1]))
		if (temp < 1):
			rho2_array[j] = temp

for i in range(len(r)):
	for j in range(M):
		x = mesh_coordinates[j][0]
		y = mesh_coordinates[j][1]
		temp = g_r(dist(x, y, r[i][0], r[i][1]))
		if (temp > 0):
			rho3_array[j] = temp

rho =  Function(VV)
rho2 = Function(V)  # Structural materials
rho3 = Function(V)  # Responsive materials

rho2.dat.data[:] = rho2_array
rho3.dat.data[:] = rho3_array

rho2 = Constant(0.0)
rho3 = Constant(0.0)

rho = as_vector([rho2, rho3])
rho = interpolate(rho, VV)

rho_initial = Function(V)
rho_initial =  rho3 - rho2
rho_initial = interpolate(rho_initial, V)
File(options.output + '/rho_initial.pvd').write(rho_initial)
# Create initial design
###### End Initial Design #####

# Define the constant parameters used in the problem
kappa = options.kappa
cw = pi/8  # Normalization parameter
lagrange_v = 1.0e-8
lagrange_r1 = options.lagrange_r1
lagrange_r2 = options.lagrange_r2

delta = 1.0e-3
epsilon = options.epsilon

kappa_d_e = kappa / (epsilon * cw)
kappa_m_e = kappa * epsilon / cw

f = Constant((0, -2))

e1 = as_vector((1, 0))
e2 = as_vector((0, 1))

# Young's modulus of the beam and poisson ratio
E_v = delta
E_r1 = options.er1modulus
E_r2 = options.er2modulus
nu = 0.3 #nu poisson ratio

mu_v = E_v/(2 * (1 + nu))
lambda_v = (E_v * nu)/((1 + nu) * (1 - 2 * nu))

mu_r1 = E_r1/(2 * (1 + nu))
lambda_r1 = (E_r1 * nu)/((1 + nu) * (1 - 2 * nu))

mu_r2 = E_r2/(2 * (1 + nu))
lambda_r2 = (E_r2 * nu)/((1 + nu) * (1 - 2 * nu))

def v_v(rho):
	return 1 - rho.sub(0) - rho.sub(1)

def v_r1(rho):
	return rho.sub(0)

def v_r2(rho):
	return rho.sub(1)

# Define h(x)=x^2
def h_v(rho):
	return pow((1 - rho.sub(0) - rho.sub(1)), options.power_p)

# Define h(x)=x^2
def h_r1(rho):
	return pow(rho.sub(0), options.power_p)

# Define h(x)=x^2
def h_r2(rho):
	return pow(rho.sub(1), options.power_p)

# Define W(x) function
def W(rho):
	return (rho.sub(0) + rho.sub(1)) * (1 - rho.sub(0)) * (1 - rho.sub(1))

# Define stress and strain tensors
def epsilon(u):
    return 0.5 * (grad(u) + grad(u).T)

# Residual strains
epsilon_star_1 =  -1 * outer(e1, e1) + outer(e2, e2)
epsilon_star_2 =  outer(e1, e1) + -1 * outer(e2, e2)

def sigma_star_1(Id):
	return lambda_r1 * tr(epsilon_star_1) * Id + 2 * mu_r1 * epsilon_star_1

def sigma_star_2(Id):
	return lambda_r2 * tr(epsilon_star_2) * Id + 2 * mu_r2 * epsilon_star_2

def sigma_v(u, Id):
    return lambda_v * tr(epsilon(u)) * Id + 2 * mu_v * epsilon(u)

def sigma_r1(u, Id):
    return lambda_r1 * tr(epsilon(u)) * Id + 2 * mu_r1 * epsilon(u)

def sigma_r2(u, Id):
    return lambda_r2 * tr(epsilon(u)) * Id + 2 * mu_r2 * epsilon(u)

# Define test function and beam displacement
v = TestFunction(VV)
u = Function(VV)
p = Function(VV)

# The left side of the beam is clamped
bcs = DirichletBC(VV, Constant((0, 0)), 7)

# Define the objective function
func1 = inner(f, u) * ds(8)
func2 = kappa_d_e * W(rho) * dx

func3_sub1 = inner(grad(v_v(rho)), grad(v_v(rho))) * dx
func3_sub2 = inner(grad(v_r1(rho)), grad(v_r1(rho))) * dx
func3_sub3 = inner(grad(v_r2(rho)), grad(v_r2(rho))) * dx

func3 = kappa_m_e * (func3_sub1 + func3_sub2 + func3_sub3)
func4 = lagrange_v * v_v(rho) * dx  # Void material
func5 = lagrange_r1 * v_r1(rho) * dx  # Structural material
func6 = lagrange_r2 * v_r2(rho) * dx  # Responsive material

J = func1 + func2 + func3 + func4 + func5 +  func6

# Define the weak form for forward PDE
a_forward_v = h_v(rho) * inner(sigma_v(u, Id), epsilon(v)) * dx
a_forward_r1 = h_r1(rho) * inner(sigma_r1(u, Id), epsilon(v)) * dx
a_forward_r2 = h_r2(rho) * inner(sigma_r2(u, Id), epsilon(v)) * dx
a_forward = a_forward_v + a_forward_r1 + a_forward_r2

L_forward_r1 = h_r1(rho) * inner(sigma_star_1(Id), epsilon(v)) * dx
L_forward_r2 = h_r2(rho) * inner(sigma_star_2(Id), epsilon(v)) * dx
L_forward = L_forward_r1 +  L_forward_r2
R_fwd = a_forward - L_forward

# Define the Lagrangian
a_lagrange_v = h_v(rho) * inner(sigma_v(u, Id), epsilon(p)) * dx
a_lagrange_r1 = h_r1(rho) * inner(sigma_r1(u, Id), epsilon(p)) * dx
a_lagrange_r2 = h_r2(rho) * inner(sigma_r2(u, Id), epsilon(p)) * dx
a_lagrange   = a_lagrange_v + a_lagrange_r1 + a_lagrange_r2

L_lagrange_r1 = h_r1(rho) * inner(sigma_star_1(Id), epsilon(p)) * dx
L_lagrange_r2 = h_r2(rho) * inner(sigma_star_2(Id), epsilon(p)) * dx
L_lagrange = L_lagrange_r1 + L_lagrange_r2
R_lagrange = a_lagrange - L_lagrange
L = J + R_lagrange


# Define the weak form for adjoint PDE
a_adjoint_v = h_v(rho) * inner(sigma_v(v, Id), epsilon(p)) * dx
a_adjoint_r1 = h_r1(rho) * inner(sigma_r1(v, Id), epsilon(p)) * dx
a_adjoint_r2 = h_r2(rho) * inner(sigma_r2(v, Id), epsilon(p)) * dx
a_adjoint = a_adjoint_v + a_adjoint_r1 + a_adjoint_r2

L_adjoint = inner(f, v) * ds(8)
R_adj = a_adjoint + L_adjoint


def FormObjectiveGradient(tao, x, G):

	i = tao.getIterationNumber()
	if (i%10) == 0:
		rho_i = Function(V)
		rho_i = rho.sub(1) - rho.sub(0)
		rho_i = interpolate(rho_i, V)
		File(options.output + '/rho-{}.pvd'.format(i)).write(rho_i)

	with rho.dat.vec as rho_vec:
		rho_vec.set(0.0)
		rho_vec.axpy(1.0, x)

	volume = assemble(rho.sub(0) * dx) * 3
	print("The volume fraction(Vr1) is {}".format(volume))

	volume = assemble(rho.sub(1) * dx) * 3
	print("The volume fraction(Vr2) is {}".format(volume))

	# Solve forward PDE
	solve(R_fwd == 0, u, bcs = bcs)

	# Solve adjoint PDE
	solve(R_adj == 0, p, bcs = bcs)

	dJdrho2 = assemble(derivative(L, rho.sub(0)))
	dJdrho3 = assemble(derivative(L, rho.sub(1)))

	dJdrho2_array = dJdrho2.vector().array()
	dJdrho3_array = dJdrho3.vector().array()

	print(dJdrho3_array)

	N = M * 2
	index_2 = []
	index_3 = []

	for i in range(N):
		if (i%2) == 0:
			index_2.append(i)
		if (i%2) == 1:
			index_3.append(i)

	G.setValues(index_2, dJdrho2_array)
	G.setValues(index_3, dJdrho3_array)

	# print(G.view())

	f_val = assemble(J)
	return f_val

# Setting lower and upper bounds
lb = as_vector((0, 0))
ub = as_vector((1, 1))
lb = interpolate(lb, VV)
ub = interpolate(ub, VV)

with lb.dat.vec as lb_vec:
	rho_lb = lb_vec

with ub.dat.vec as ub_vec:
	rho_ub = ub_vec

# Setting TAO solver
tao = PETSc.TAO().create(PETSc.COMM_SELF)
# tao.setType('bncg')
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

rho_final = Function(V)
rho_final = rho.sub(1) - rho.sub(0)
rho_final = interpolate(rho_final, V)
File(options.output + '/rho-final.pvd').write(rho_final)
File(options.output + '/rho-final-rho2.pvd').write(rho.sub(0))
File(options.output + '/rho-final-rho3.pvd').write(rho.sub(1))
File(options.output + '/u.pvd').write(u)

end = time.time()
print("\nExecution time (in seconds):", (end - start))
