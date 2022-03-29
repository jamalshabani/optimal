from firedrake import *
from petsc4py import PETSc
import numpy as np

# Import gmesh
mesh = Mesh("main.msh")
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
	if (s < 0.02):
		return 1
	if (0.02 < s < 0.03):
		return (cos((s - 0.02)*pi/0.01) + 1)/2
	if (s > 0.03):
		return 0

V = FunctionSpace(mesh, 'CG', 1)
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

rho = interpolate(rho, V)

File("initial/rho.pvd").write(rho)
