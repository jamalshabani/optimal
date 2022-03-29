#!/usr/bin/env python3
# Blocking load implementation
def parse():
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('-tao_monitor', '--tao_monitor', action='store_true', help = 'TAO monitor')
	parser.add_argument('-tao_max_it', '--tao_max_it', type = int, default = 100, help = 'Number of TAO iterations')
	parser.add_argument('-l', '--lagrange', type = float, default = 1.0, help = 'Lagrange multiplier')
	parser.add_argument('-k', '--kappa', type = float, default = 1.0e-2, help = 'Weight of Modica-Mortola')
	parser.add_argument('-e', '--epsilon', type = float, default = 5.0e-3, help = 'Phase-field regularization parameter')
	parser.add_argument('-o', '--output', type = str, default = 'output1', help = 'Output folder')
	parser.add_argument('-v', '--volume_fraction', type = float, default = 0.25, help = 'Volume fraction for responsive material')
	parser.add_argument('-m', '--mesh', type = str, default = 'main.msh', help = 'Dimensions of meshed beam')
	parser.add_argument('-er', '--ermodulus', type = float, default = 0.1, help = 'Elastic Modulus for responsive material')
	parser.add_argument('-es', '--esmodulus', type = float, default = 1.0, help = 'Elastic Modulus for structural material')
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
VV = VectorFunctionSpace(mesh, 'CG', 1)

# Initial Design
mesh_coordinates = mesh.coordinates.dat.data[:]

# Define the points
r1 = [[0.1, 0], [0.3, 0], [0.5, 0], [0.7, 0], [0.9, 0]]

r2 = [[0, 1/12], [0.2, 1/12], [0.4, 1/12], [0.6, 1/12], [0.8, 1/12], [1.0, 1/12]]

r3 = [[0.1, 2/12], [0.3, 2/12], [0.5, 2/12], [0.7, 2/12], [0.9, 2/12]]

r4 = [[0, 3/12], [0.2, 3/12], [0.4, 3/12], [0.6, 3/12], [0.8, 3/12], [1.0, 3/12]]

r5 = [[0.1, 4/12], [0.3, 4/12], [0.5, 4/12], [0.7, 4/12], [0.9, 4/12]]

r = r1 + r2 + r3 + r4 + r5

def dist(x, y, a, b):
	return sqrt((x - a)**2 + (y - b)**2)

def g(s):
	if (s < 0.04):
		return 1
	if (0.04 < s < 0.05):
		return (cos((s - 0.04)*pi/0.01) + 1)/2
	if (s > 0.05):
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
rho.dat.data[:] = rho_array
rho = Constant(0.4)
rho = interpolate(rho, V)
File(options.output + '/rho_initial.pvd').write(rho)
###### End Initial Design #####

# Define the constant parameters used in the problem
kappa = options.kappa
cw = pi/8  # Normalization parameter
lagrange = options.lagrange # Lagrange multiplier for Volume constraint
epsilon = options.epsilon
volume_fraction_0 = options.volume_fraction
volume_tol = 1e-3

kappa_d_e = kappa / (epsilon * cw)
kappa_m_e = kappa * epsilon / cw

f = Constant((0, -2))

e1 = as_vector((1, 0))
e2 = as_vector((0, 1))

# Young's modulus of the beam and poisson ratio
E_s = options.esmodulus
E_r = options.ermodulus
nu = 0.3 #nu poisson ratio

mu_s = E_s/(2 * (1 + nu))
lambda_s = (E_s * nu)/((1 + nu) * (1 - 2 * nu))

mu_r = E_r/(2 * (1 + nu))
lambda_r = (E_r * nu)/((1 + nu) * (1 - 2 * nu))


def v_s(rho):
	return (1 - rho)

def v_r(rho):
	return rho

# Define h(x)=x^p where p >=2
def h_s(rho):
	return pow((1 - rho), options.power_p)

# Define h(x)=x^p where p >=2
def h_r(rho):
	return pow(rho, options.power_p)

# Define W(x) function
def W(rho):
	return rho * (1 - rho)

# Define stress and strain tensors
def epsilon(u):
	return 0.5 * (grad(u) + grad(u).T)

# Residual strain
epsilon_star =  -1 * outer(e1, e1) + outer(e2, e2)

def sigma_star(Id):
	return lambda_r * tr(epsilon_star) * Id + 2 * mu_r * epsilon_star

def sigma_s(u, Id):
	return lambda_s * tr(epsilon(u)) * Id + 2 * mu_s * epsilon(u)

def sigma_r(u, Id):
	return lambda_r * tr(epsilon(u)) * Id + 2 * mu_r * epsilon(u)

# Define test function and beam displacement
v = TestFunction(VV)
u = Function(VV)
p = Function(VV)

# The left side of the beam is clamped
bcs = DirichletBC(VV, Constant((0, 0)), 7)

# Define the objective function
func1 = inner(f, u) * ds(8)
func2 = kappa_d_e * W(rho) * dx

func3_sub1 = inner(grad(v_s(rho)), grad(v_s(rho))) * dx
func3_sub2 = inner(grad(v_r(rho)), grad(v_r(rho))) * dx

func3 = kappa_m_e * (func3_sub1 + func3_sub2)
func4 = lagrange * v_r(rho) * dx

J = func1 + func2 + func3 + func4

# Define the weak form for forward PDE
a_forward_s = h_s(rho) * inner(sigma_s(u, Id), epsilon(v)) * dx
a_forward_r = h_r(rho) * inner(sigma_r(u, Id), epsilon(v)) * dx
a_forward = a_forward_s + a_forward_r

L_forward = h_r(rho) * inner(sigma_star(Id), epsilon(v)) * dx
R_fwd = a_forward - L_forward

# Define the Lagrangian
a_lagrange_s = h_s(rho) * inner(sigma_s(u, Id), epsilon(p)) * dx
a_lagrange_r = h_r(rho) * inner(sigma_r(u, Id), epsilon(p)) * dx
a_lagrange   = a_lagrange_s + a_lagrange_r

L_lagrange = h_r(rho) * inner(sigma_star(Id), epsilon(p)) * dx
R_lagrange = a_lagrange - L_lagrange
L = J + R_lagrange


# Define the weak form for adjoint PDE
a_adjoint_s = h_s(rho) * inner(sigma_s(v, Id), epsilon(p)) * dx
a_adjoint_r = h_r(rho) * inner(sigma_r(v, Id), epsilon(p)) * dx
a_adjoint = a_adjoint_s + a_adjoint_r

L_adjoint = inner(f, v) * ds(8)
R_adj = a_adjoint + L_adjoint


def FormObjectiveGradient(tao, x, G):

	i = tao.getIterationNumber()
	if (i%10) == 0:
		File(options.output + '/rho-{}.pvd'.format(i)).write(rho)

	if (options.mesh == '1_to_3_mesh.msh'):
		volume_fraction = assemble(rho * dx) * 3
	elif (options.mesh == '1_to_1_mesh.msh'):
		volume_fraction = assemble(rho * dx)
	else:
		volume_fraction = assemble(rho * dx) * 6

	print("The volume fraction(Vr) is {}".format(volume_fraction))

	# if (abs(volume_fraction - volume_fraction_0) < volume_tol):
	# 	tao.setConvergedReason(8)

	with rho.dat.vec as rho_vec:
		rho_vec.set(0.0)
		rho_vec.axpy(1.0, x)

	# Solve forward PDE
	solve(R_fwd == 0, u, bcs = bcs)

	# Solve adjoint PDE
	solve(R_adj == 0, p, bcs = bcs)

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

File(options.output + '/final-rho.pvd').write(rho)
File(options.output + '/displacement.pvd').write(u)

end = time.time()
print("\nExecution time (in seconds):", (end - start))
