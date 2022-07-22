def parse():
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('-tao_type', '--tao_type', type = str, default = 'blmvm', help = 'TAO algorithm type')
	parser.add_argument('-tao_monitor', '--tao_monitor', action='store_true', help = 'TAO monitor')
	parser.add_argument('-tao_max_it', '--tao_max_it', type = int, default = 100, help = 'Number of TAO iterations')
	parser.add_argument('-tao_gatol', '--tao_gatol', type = float, default = 1.0e-7, help = 'Stop if norm of gradient is less than this')
	parser.add_argument('-tao_grtol', '--tao_grtol', type = float, default = 1.0e-7, help = 'Stop if relative norm of gradient is less than this')
	parser.add_argument('-tao_gttol', '--tao_gttol', type = float, default = 1.0e-7, help = 'Stop if norm of gradient is reduced by this factor')
	parser.add_argument('-ls', '--lagrange_s', type = float, default = 1.0, help = 'Lagrange multiplier for structural material')
	parser.add_argument('-lr', '--lagrange_r', type = float, default = 0.1, help = 'Lagrange multiplier for responsive material')
	parser.add_argument('-vs', '--volume_s', type = float, default = 0.4, help = 'Volume percentage for structural material')
	parser.add_argument('-vr', '--volume_r', type = float, default = 0.3, help = 'Volume percentage for responsive material')
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
VV = VectorFunctionSpace(mesh, 'CG', 1, dim = 2)

# Create initial design
###### Begin Initial Design #####
mesh_coordinates = mesh.coordinates.dat.data[:]

M = len(mesh_coordinates)
rho2_array = np.ones(M) # Blue material
rho3_array = np.zeros(M) # Red material

rho =  Function(VV, name = "Design variable")
rho2 = Function(V, name = "Structural material")  # Structural material 1(Blue)
rho3 = Function(V, name = "Responsive material")  # Responsive material 2(Red)

x, y = SpatialCoordinate(mesh)
rho2 = Constant(0.4)
rho3 = Constant(0.4)
# rho2 = 0.75 + 0.75 * sin(4*pi*x) * sin(8*pi*y)
# rho3 = 0.50 + 0.50 * sin(4*pi*x) * sin(8*pi*y)

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
lagrange_s = options.lagrange_s
lagrange_r = options.lagrange_r
volume_s = options.volume_s
volume_r = options.lagrange_r

# Total volume of the domain |omega|
omega = assemble(interpolate(Constant(1.0), V) * dx)

delta = 1.0e-3
epsilon = options.epsilon
kappa_d_e = kappa / (epsilon * cw)
kappa_m_e = kappa * epsilon / cw

f = Constant((0, 0.0))
u_star = Constant((0, -1.0))

# Young's modulus of the beam and poisson ratio
E_v = delta
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

# Residual strain
e1 = as_vector((1, 0))
epsilon_star =  0.5 * (Id - 2 * outer(e1, e1))

def sigma_star(Id):
	return lambda_r * tr(epsilon_star) * Id + 2 * mu_r * epsilon_star

def sigma_v(u, Id):
    return lambda_v * tr(epsilon(u)) * Id + 2 * mu_v * epsilon(u)

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

func2_sub1 = inner(grad(v_v(rho)), grad(v_v(rho))) * dx
func2_sub2 = inner(grad(v_s(rho)), grad(v_s(rho))) * dx
func2_sub3 = inner(grad(v_r(rho)), grad(v_r(rho))) * dx

func2 = kappa_m_e * (func2_sub1 + func2_sub2 + func2_sub3)
func3 = lagrange_s * (v_s(rho) - volume_s * omega) * dx  # Responsive material 1(Blue)
func4 = lagrange_r * (v_r(rho) - volume_r * omega) * dx  # Responsive material 2(Red)

# Objective function + Modica-Mortola functional + Volume constraint
P = func1 + func2 + func3 + func4
JJ = J + P

# Define the weak form for forward PDE
a_forward_v = h_v(rho) * inner(sigma_v(u, Id), epsilon(v)) * dx
a_forward_s = h_s(rho) * inner(sigma_s(u, Id), epsilon(v)) * dx
a_forward_r = h_r(rho) * inner(sigma_r(u, Id), epsilon(v)) * dx
a_forward = a_forward_v + a_forward_s + a_forward_r

L_forward = h_r(rho) * inner(sigma_star(Id), epsilon(v)) * dx
R_fwd = a_forward - L_forward

# Define the Lagrangian
a_lagrange_v = h_v(rho) * inner(sigma_v(u, Id), epsilon(p)) * dx
a_lagrange_s = h_s(rho) * inner(sigma_s(u, Id), epsilon(p)) * dx
a_lagrange_r = h_r(rho) * inner(sigma_r(u, Id), epsilon(p)) * dx
a_lagrange   = a_lagrange_v + a_lagrange_s + a_lagrange_r

L_lagrange = h_r(rho) * inner(sigma_star(Id), epsilon(p)) * dx
R_lagrange = a_lagrange - L_lagrange
L = JJ - R_lagrange


# Define the weak form for adjoint PDE
a_adjoint_v = h_v(rho) * inner(sigma_v(v, Id), epsilon(p)) * dx
a_adjoint_s = h_s(rho) * inner(sigma_s(v, Id), epsilon(p)) * dx
a_adjoint_r = h_r(rho) * inner(sigma_r(v, Id), epsilon(p)) * dx
a_adjoint = a_adjoint_v + a_adjoint_s + a_adjoint_r

L_adjoint = inner(u - u_star, v) * dx(4)
R_adj = a_adjoint - L_adjoint


def FormObjectiveGradient(tao, x, G):

	i = tao.getIterationNumber()
	if (i%10) == 0:
		rho_i = Function(V)
		rho_i = rho.sub(1) - rho.sub(0)
		rho_i = interpolate(rho_i, V)
		# File(options.output + '/rho-{}.pvd'.format(i)).write(rho_i)
		# File(options.output + '/u-{}.pvd'.format(i)).write(u)
		File(options.output + '/beam-{}.pvd'.format(i)).write(rho_i, u)


	with rho.dat.vec as rho_vec:
		rho_vec.set(0.0)
		rho_vec.axpy(1.0, x)

	# Solve forward PDE
	solve(R_fwd == 0, u, bcs = bcs, solver_parameters={'snes_max_it': 500})

	# Solve adjoint PDE
	solve(R_adj == 0, p, bcs = bcs, solver_parameters={'snes_max_it': 500})

	objective_value = assemble(J)
	print("The value of objective function is {}".format(objective_value))

	volume_s = assemble(v_s(rho) * dx)/omega
	print("The volume fraction(Vs) is {}".format(volume_s))

	volume_r = assemble(v_r(rho) * dx)/omega
	print("The volume fraction(Vr) is {}".format(volume_r))
	print(" ")

	dJdrho2 = assemble(derivative(L, rho.sub(0)))
	dJdrho3 = assemble(derivative(L, rho.sub(1)))

	dJdrho2_array = dJdrho2.vector().array()
	dJdrho3_array = dJdrho3.vector().array()

	# print(dJdrho2_array)

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
	f_val = assemble(JJ)
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
rho_final = rho.sub(1) - rho.sub(0)
rho_final = interpolate(rho_final, V)
File(options.output + '/rho-final.pvd').write(rho_final)
File(options.output + '/rho-final-rho2.pvd').write(rho.sub(0))
File(options.output + '/rho-final-rho3.pvd').write(rho.sub(1))
File(options.output + '/displacement.pvd').write(u)
File(options.output + '/beam-final.pvd').write(rho_final, u)

end = time.time()
print("\nExecution time (in seconds):", (end - start))
