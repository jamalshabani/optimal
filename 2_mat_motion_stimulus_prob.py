#!/usr/bin/env python3
def parse():
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('-tao_type', '--tao_type', type = str, default = 'blmvm', help = 'TAO algorithm type')
	parser.add_argument('-tao_monitor', '--tao_monitor', action='store_true', help = 'TAO monitor')
	parser.add_argument('-tao_max_it', '--tao_max_it', type = int, default = 100, help = 'Number of TAO iterations')
	parser.add_argument('-tao_gatol', '--tao_gatol', type = float, default = 1.0e-7, help = 'Stop if norm of gradient is less than this')
	parser.add_argument('-tao_grtol', '--tao_grtol', type = float, default = 1.0e-7, help = 'Stop if relative norm of gradient is less than this')
	parser.add_argument('-tao_gttol', '--tao_gttol', type = float, default = 1.0e-7, help = 'Stop if norm of gradient is reduced by this factor')
	parser.add_argument('-l', '--lagrange', type = float, default = 0.1, help = 'Lagrange multiplier for responsive material')
	parser.add_argument('-v', '--volume', type = float, default = 0.4, help = 'Volume percentage for structural material')
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
mesh = Mesh(options.mesh)
Id = Identity(mesh.geometric_dimension()) #Identity tensor

# Define the function spaces
V = FunctionSpace(mesh, 'CG', 1)
VV = VectorFunctionSpace(mesh, 'CG', 1)

mesh_coordinates = mesh.coordinates.dat.data[:]

M = len(mesh_coordinates)
print(M)

# Initial Design and stimulus
rho = Function(V)
# rho = 0.75 + 0.75 * sin(4*pi*x) * sin(8*pi*y)
s = 1.0
rho = Constant(0.4)
rho = interpolate(rho, V)
rhos = as_vector([rho, s])
rhos = interpolate(rhos, VV)
print(type(rhos))
print(type(rhos.sub(0)))
print(type(rhos.sub(1)))
File(options.output + '/rho_initial.pvd').write(rho)
###### End Initial Design #####

# Define the constant parameters used in the problem
kappa = options.kappa
cw = pi/8  # Normalization parameter
lagrange = options.lagrange # Lagrange multiplier for Volume constraint
epsilon = options.epsilon
volume = options.volume

# Total volume of the domain |omega|
omega = assemble(interpolate(Constant(1.0), V) * dx)

kappa_d_e = kappa / (epsilon * cw)
kappa_m_e = kappa * epsilon / cw

f = Constant((0, -1))
u_star = Constant((0, 1))

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


def v_s(rhos):
	return (1 - rho.sub(0))

def v_r(rhos):
	return rho.sub(0)

def a_s(rhos):
	return 2 * pow(rhos.sub(1), 2) - 1

def b_s(rhos):
	return 2 * rhos.sub(1) * sqrt((1 - rhos.sub(1)))

# Define h(x)=x^p where p >=2
def h_s(rhos):
	return pow((1 - rho.sub(0)), options.power_p)

# Define h(x)=x^p where p >=2
def h_r(rhos):
	return pow(rho.sub(0), options.power_p)

# Define W(x) function
def W(rhos):
	return rho.sub(0) * (1 - rho.sub(0))

# Define stress and strain tensors
def epsilon(u):
	return 0.5 * (grad(u) + grad(u).T)

# Residual strain
A = outer(e1, e1) - outer(e2, e2)
B = outer(e1, e2) + outer(e2, e1)

def sigma_A(A, Id):
	return lambda_r * tr(A) * Id + 2 * mu_r * A

def sigma_B(B, Id):
	return lambda_r * tr(B) * Id + 2 * mu_r * B

def sigma_s(u, Id):
	return lambda_s * tr(epsilon(u)) * Id + 2 * mu_s * epsilon(u)

def sigma_r(u, Id):
	return lambda_r * tr(epsilon(u)) * Id + 2 * mu_r * epsilon(u)

# Define test function and beam displacement
v = TestFunction(VV)
u = Function(VV, name = "Displacement")
p = Function(VV, name = "Adjoint variable")

# The left side of the beam is clamped
bcs = DirichletBC(VV, Constant((0, 0)), 7)

# Define the objective function
J = 0.5 * inner(u - u_star, u - u_star) * dx(4)
func1 = kappa_d_e * W(rhos) * dx

func2_sub1 = inner(grad(v_s(rhos)), grad(v_s(rhos))) * dx
func2_sub2 = inner(grad(v_r(rhos)), grad(v_r(rhos))) * dx

func2 = kappa_m_e * (func2_sub1 + func2_sub2)
func3 = lagrange * (v_r(rhos) - volume * omega) * dx

JJ = J + func1 + func2 + func3

# Define the weak form for forward PDE
a_forward_s = h_s(rhos) * inner(sigma_s(u, Id), epsilon(v)) * dx
a_forward_r = h_r(rhos) * inner(sigma_r(u, Id), epsilon(v)) * dx
a_forward_A = a_s(rhos) * h_r(rhos) * inner(sigma_A(A, Id), epsilon(v)) * dx
a_forward_B = b_s(rhos) * h_r(rhos) * inner(sigma_B(B, Id), epsilon(v)) * dx
a_forward = a_forward_s + a_forward_r + a_forward_A + a_forward_B

L_forward = inner(f, v) * ds(8)
R_fwd = a_forward - L_forward

# Define the Lagrangian
a_lagrange_s = h_s(rhos) * inner(sigma_s(u, Id), epsilon(p)) * dx
a_lagrange_r = h_r(rhos) * inner(sigma_r(u, Id), epsilon(p)) * dx
a_lagrange_A = a_s(rhos) * h_r(rhos) * inner(sigma_A(A, Id), epsilon(p)) * dx
a_lagrange_B = b_s(rhos) * h_r(rhos) * inner(sigma_B(B, Id), epsilon(p)) * dx
a_lagrange   = a_lagrange_s + a_lagrange_r + a_lagrange_A + a_lagrange_B

L_lagrange = inner(f, p) * ds(8)
R_lagrange = a_lagrange - L_lagrange
L = JJ + R_lagrange


# Define the weak form for adjoint PDE
a_adjoint_s = h_s(rhos) * inner(sigma_s(v, Id), epsilon(p)) * dx
a_adjoint_r = h_r(rhos) * inner(sigma_r(v, Id), epsilon(p)) * dx
a_adjoint = a_adjoint_s + a_adjoint_r

L_adjoint = inner(u - u_star, v) * dx(4)
R_adj = a_adjoint + L_adjoint


def FormObjectiveGradient(tao, x, G):

	i = tao.getIterationNumber()
	if (i%10) == 0:
		File(options.output + '/beam-{}.pvd'.format(i)).write(rhos.sub(0), u)

	objective_value = assemble(J)
	print("The value of objective function is {}".format(objective_value))

	volume_fraction = assemble(rhos.sub(0) * dx) * 3
	print("The volume fraction(Vr) is {}".format(volume_fraction))
	print(" ")

	with rhos.dat.vec as rho_vec:
		rho_vec.set(0.0)
		rho_vec.axpy(1.0, x)

	# Solve forward PDE
	solve(R_fwd == 0, u, bcs = bcs)

	# Solve adjoint PDE
	solve(R_adj == 0, p, bcs = bcs)

	dJdrho = assemble(derivative(L, rhos.sub(0)))
	dJds = assemble(derivative(L, rhos.sub(1)))

	dJdrho_array = dJdrho.vector().array()
	dJds_array = dJds.vector().array()

	# print(dJdrho2_array)

	N = M * 2
	index_rho = []
	index_s = []

	for i in range(N):
		if (i%2) == 0:
			index_rho.append(i)
		if (i%2) == 1:
			index_s.append(i)

	G.setValues(index_rho, dJdrho_array)
	G.setValues(index_s, dJds_array)

	# dJdrho = assemble(derivative(L, rho))
	# with dJdrho.dat.vec as dJdrho_vec:
	# 	G.set(0.0)
	# 	G.axpy(1.0, dJdrho_vec)

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
tao.setType('blmvm')
tao.setObjectiveGradient(FormObjectiveGradient, None)
tao.setVariableBounds(rho_lb, rho_ub)
tao.setFromOptions()

# Initial design guess
with rhos.dat.vec as rhos_vec:
	x = rhos_vec.copy()

# Solve the optimization problem
tao.solve(x)
tao.destroy()

# Recover the final solution
with rhos.dat.vec as rhos_vec:
	rhos_vec = x.copy()

# Saving files for viewing with Paraview
rho_final = Function(V, name = "Design variable")
rho_final = rhos.sub(0)
s_final = rhos.sub(1)
rho_final = interpolate(rho_final, V)
File(options.output + '/rho-final.pvd').write(rho_final)
File(options.output + '/rho-final-rho.pvd').write(rho.sub(0))
File(options.output + '/s-final.pvd').write(rho.sub(1))
File(options.output + '/displacement.pvd').write(u)
File(options.output + '/beam-final.pvd').write(rho_final, u)

end = time.time()
print("\nExecution time (in seconds):", (end - start))
