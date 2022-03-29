from firedrake import *
from petsc4py import PETSc
import numpy as np

# Import gmesh
mesh = Mesh("main.msh")
mesh_coordinates = mesh.coordinates.dat.data[:]

def g(s):
    if (s >= 1/6):
        return 1
    if (s < 1/6):
        return 0

V = FunctionSpace(mesh, 'CG', 1)
M = len(mesh_coordinates)

rho_array = np.zeros(M)

for j in range(M):
	y = mesh_coordinates[j][1]
	temp = g(y)
	if (temp > 0):
		rho_array[j] = temp


rho = Function(V)
rho.dat.data[:] = rho_array

rho = interpolate(rho, V)

File("initial/rho.pvd").write(rho)
