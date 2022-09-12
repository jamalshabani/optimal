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
	parser.add_argument('-m', '--mesh', type = str, default = 'motion_mesh.msh', help = 'Dimensions of meshed beam')
	parser.add_argument('-es', '--esmodulus', type = float, default = 0.1, help = 'Elastic Modulus for structural material')
	parser.add_argument('-er', '--ermodulus', type = float, default = 1.0, help = 'Elastic Modulus for responsive material')
	parser.add_argument('-p', '--power_p', type = float, default = 2.0, help = 'Power for elasticity interpolation')
	options = parser.parse_args()
	return options

options = parse()

from firedrake import *
from petsc4py import PETSc
import time

start = time.time()

# Import gmesh
mesh = Mesh(options.mesh)
Id = Identity(mesh.geometric_dimension()) #Identity tensor

# Define the function spaces
V = FunctionSpace(mesh, 'CG', 1)
VV = VectorFunctionSpace(mesh, 'CG', 1, dim = 2)

# Initial Design and stimulus
rho = Function(V)
x, y = SpatialCoordinate(mesh)
# rho = 0.75 + 0.75 * sin(4*pi*x) * sin(8*pi*y)
rho = Constant(0.4)
rho = interpolate(rho, V)
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

f = Constant((0, -1.0))
u_star = Constant((0, 1.0))

# Young's modulus of the beam and poisson ratio
E_s = options.esmodulus
E_r = options.ermodulus
nu = 0.3 # Poisson ratio

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
s = Constant(1.0)
e1 = as_vector((1, 0)) # Direction of responsive material
epsilon_star =  Id - 2 * outer(e1, e1)

def sigma_a(A, Id):
	return lambda_r * tr(A) * Id + 2 * mu_r * A

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
func1 = kappa_d_e * W(rho) * dx

func2_sub1 = inner(grad(v_s(rho)), grad(v_s(rho))) * dx
func2_sub2 = inner(grad(v_r(rho)), grad(v_r(rho))) * dx

func2 = kappa_m_e * (func2_sub1 + func2_sub2)
func3 = lagrange * (v_r(rho) - volume * omega) * dx

JJ = J + func1 + func2 + func3

# Define the weak form for forward PDE
a_forward_s = h_s(rho) * inner(sigma_s(u, Id), epsilon(v)) * dx
a_forward_r = h_r(rho) * inner(sigma_r(u, Id), epsilon(v)) * dx
a_forward = a_forward_s + a_forward_r

L_forward = inner(f, v) * ds(8) + h_r(rho) * inner(sigma_a(epsilon_star, Id), epsilon(v)) * dx
R_fwd = a_forward - L_forward

# Define the Lagrangian
a_lagrange_s = h_s(rho) * inner(sigma_s(u, Id), epsilon(p)) * dx
a_lagrange_r = h_r(rho) * inner(sigma_r(u, Id), epsilon(p)) * dx
a_lagrange   = a_lagrange_s + a_lagrange_r

L_lagrange = inner(f, p) * ds(8) + h_r(rho) * inner(sigma_a(epsilon_star, Id), epsilon(p)) * dx
R_lagrange = a_lagrange - L_lagrange
L = JJ - R_lagrange


# Define the weak form for adjoint PDE
a_adjoint_s = h_s(rho) * inner(sigma_s(v, Id), epsilon(p)) * dx
a_adjoint_r = h_r(rho) * inner(sigma_r(v, Id), epsilon(p)) * dx
a_adjoint = a_adjoint_s + a_adjoint_r

L_adjoint = inner(u - u_star, v) * dx(4)
R_adj = a_adjoint - L_adjoint


def FormObjectiveGradient(tao, x, G):

	i = tao.getIterationNumber()
	if (i%10) == 0:
		File(options.output + '/beam-{}.pvd'.format(i)).write(rho, u)

	volume_fraction = assemble(rho * dx) * 3
	print("\nThe volume fraction(Vr) is {:0.6f}".format(volume_fraction))

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

	f_val = assemble(JJ)
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

# Saving files for viewing with Paraview
rho_final = Function(V, name = "Design variable")
File(options.output + '/rho-final.pvd').write(rho)
File(options.output + '/displacement.pvd').write(u)
File(options.output + '/beam-final.pvd').write(rho, u)
print(" ")
print("-------------------------------------------------------")

objective_value = assemble(J)
print("The final value of objective function is {:0.6f}".format(objective_value))

end = time.time()
print("\nExecution time (in seconds): {:0.6f}".format((end - start)))
print(" ")
