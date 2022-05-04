def parse():
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('-tao_monitor', '--tao_monitor', action='store_true', help = 'TAO monitor')
	parser.add_argument('-tao_max_it', '--tao_max_it', type = int, default = 100, help = 'Number of TAO iterations')
	parser.add_argument('-ls', '--lagrange_s', type = float, default = 5.0, help = 'Lagrange multiplier for structural material')
	parser.add_argument('-lr', '--lagrange_r', type = float, default = 0.5, help = 'Lagrange multiplier for responsive material')
	parser.add_argument('-k', '--kappa', type = float, default = 1.0e-2, help = 'Weight of Modica-Mortola')
	parser.add_argument('-e', '--epsilon', type = float, default = 5.0e-3, help = 'Phase-field regularization parameter')
	parser.add_argument('-o', '--output', type = str, default = 'output1', help = 'Output folder')
	parser.add_argument('-m', '--mesh', type = str, default = 'main.msh', help = 'Dimensions of meshed beam')
	parser.add_argument('-es', '--esmodulus', type = float, default = 0.1, help = 'Elastic Modulus for structural material')
	parser.add_argument('-er', '--ermodulus', type = float, default = 1.0, help = 'Elastic Modulus for responsive material')
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
mesh = Mesh("main.msh")
Id = Identity(mesh.geometric_dimension()) #Identity tensor

# Define the function spaces
V = FunctionSpace(mesh, 'CG', 1)
VV = VectorFunctionSpace(mesh, 'CG', 1, dim = 2)

# Create initial design
###### Begin Initial Design #####
mesh_coordinates = mesh.coordinates.dat.data[:]

# Define the points
v1 = [[0.1, 1/12], [0.3, 1/12], [0.5, 1/12], [0.7, 1/12], [0.9, 1/12]]
r1 = [[0.2, 1/12], [0.4, 1/12], [0.6, 1/12], [0.8, 1/12]]

v2 = [[0.2, 2/12], [0.4, 2/12], [0.6, 2/12], [0.8, 2/12]]
r2 = [[0.1, 2/12], [0.3, 2/12], [0.5, 2/12], [0.7, 2/12], [0.9, 2/12]]

v3 = [[0.1, 3/12], [0.3, 3/12], [0.5, 3/12], [0.7, 3/12], [0.9, 3/12]]
r3 = [[0.2, 3/12], [0.4, 3/12], [0.6, 3/12], [0.8, 3/12]]

v = v1 + v2 + v3
r = r1 + r2 + r3
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

# rho2 = Constant(0.4)
# rho3 = Constant(0.4)

rho = as_vector([rho2, rho3])
rho = interpolate(rho, VV)

rho_initial = Function(V)
rho_initial =  rho3 - rho2
rho_initial = interpolate(rho_initial, V)
File("output4/rho_initial.pvd").write(rho_initial)
# Create initial design
###### End Initial Design #####

# Define the constant parameters used in the problem
kappa = options.kappa
lagrange_v = 1.0e-8
lagrange_s = options.lagrange_s
lagrange_r = options.lagrange_r
epsilon = options.epsilon

kappa_d_e = kappa / epsilon
kappa_m_e = kappa * epsilon

f = Constant((0, -1))

# Young's modulus of the beam and poisson ratio
E_v = 1.0e-3
E_s = options.esmodulus
E_r = options.ermodulus
nu = 0.3 #nu poisson ratio

mu_v = E_v/(2 * (1 + nu))
lambda_v = (E_v * nu)/((1 + nu) * (1 - 2 * nu))

mu_s = E_s/(2 * (1 + nu))
lambda_s = (E_s * nu)/((1 + nu) * (1 - 2 * nu))

mu_r = E_r/(2 * (1 + nu))
lambda_r = (E_r * nu)/((1 + nu) * (1 - 2 * nu))

def v_v(rho):
	return 1 - rho.sub(0) - rho.sub(1)

def v_s(rho):
	return rho.sub(0)

def v_r(rho):
	return rho.sub(1)

# Define h(x)=x^2
def h_v(rho):
	return pow((1 - rho.sub(0) - rho.sub(1)), options.power_p)

# Define h(x)=x^2
def h_s(rho):
	return pow(rho.sub(0), options.power_p)

# Define h(x)=x^2
def h_r(rho):
	return pow(rho.sub(1), options.power_p)

# Define W(x) function
def W(rho):
	return (rho.sub(0) + rho.sub(1)) * (1 - rho.sub(0)) * (1 - rho.sub(1))

# Define stress and strain tensors
def epsilon(u):
    return 0.5 * (grad(u) + grad(u).T)

def sigma_v(u, Id):
    return lambda_v * tr(epsilon(u)) * Id + 2 * mu_v * epsilon(u)

def sigma_s(u, Id):
    return lambda_s * tr(epsilon(u)) * Id + 2 * mu_s * epsilon(u)

def sigma_r(u, Id):
    return lambda_r * tr(epsilon(u)) * Id + 2 * mu_r * epsilon(u)

# Define test function and beam displacement
v = TestFunction(VV)
u = Function(VV)

# The left side of the beam is clamped
bcs = DirichletBC(VV, Constant((0, 0)), 7)

# Define the objective function
func1 = inner(f, u) * ds(8)
func2 = kappa_d_e * W(rho) * dx

func3_sub1 = inner(grad(v_v(rho)), grad(v_v(rho))) * dx
func3_sub2 = inner(grad(v_s(rho)), grad(v_s(rho))) * dx
func3_sub3 = inner(grad(v_r(rho)), grad(v_r(rho))) * dx

func3 = kappa_m_e * 0.5 * (func3_sub1 + func3_sub2 + func3_sub3)
func4 = lagrange_v * v_v(rho) * dx  # Void material
func5 = lagrange_s * v_s(rho) * dx  # Structural material
func6 = lagrange_r * v_r(rho) * dx  # Responsive material

J = func1 + func2 + func3 + func4 + func5 + func6

# Define the weak form for forward PDE
a_forward_v = h_v(rho) * inner(sigma_v(u, Id), epsilon(v)) * dx
a_forward_s = h_s(rho) * inner(sigma_s(u, Id), epsilon(v)) * dx
a_forward_r = h_r(rho) * inner(sigma_r(u, Id), epsilon(v)) * dx
a_forward = a_forward_v + a_forward_s + a_forward_r

L_forward = inner(f, v) * ds(8)
R_fwd = a_forward - L_forward

# Define the Lagrangian
a_lagrange_v = h_v(rho) * inner(sigma_v(u, Id), epsilon(u)) * dx
a_lagrange_s = h_s(rho) * inner(sigma_s(u, Id), epsilon(u)) * dx
a_lagrange_r = h_r(rho) * inner(sigma_r(u, Id), epsilon(u)) * dx
a_lagrange   = a_lagrange_v + a_lagrange_s + a_lagrange_r

L_lagrange = inner(f, u) * ds(8)
R_lagrange = a_lagrange - L_lagrange
L = J - R_lagrange


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
	print("The volume fraction(Vs) is {}".format(volume))

	volume = assemble(rho.sub(1) * dx) * 3
	print("The volume fraction(Vr) is {}".format(volume))

	# Solve forward PDE
	solve(R_fwd == 0, u, bcs = bcs)

	dJdrho2 = assemble(derivative(L, rho.sub(0)))
	dJdrho3 = assemble(derivative(L, rho.sub(1)))

	dJdrho2_array = dJdrho2.vector().array()
	dJdrho3_array = dJdrho3.vector().array()

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
tao.setType('bncg')
# tao.setType('blmvm')
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
