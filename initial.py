from firedrake import *
from petsc4py import PETSc
import numpy as np

# Import gmesh
mesh = Mesh("main.msh")
mesh_coordinates = mesh.coordinates.dat.data[:]


# Define the points
v1 = [[0, 0], [0.2, 0], [0.4, 0], [0.6, 0], [0.8, 0], [1.0, 0]]
r1 = [[0.1, 0], [0.3, 0], [0.5, 0], [0.7, 0], [0.9, 0]]

v2 = [[0.1, 1/12], [0.3, 1/12], [0.5, 1/12], [0.7, 1/12], [0.9, 1/12]]
r2 = [[0, 1/12], [0.2, 1/12], [0.4, 1/12], [0.6, 1/12], [0.8, 1/12], [1.0, 1/12]]

v3 = [[0, 2/12], [0.2, 2/12], [0.4, 2/12], [0.6, 2/12], [0.8, 2/12], [1.0, 2/12]]
r3 = [[0.1, 2/12], [0.3, 2/12], [0.5, 2/12], [0.7, 2/12], [0.9, 2/12]]

v4 = [[0.1, 3/12], [0.3, 3/12], [0.5, 3/12], [0.7, 3/12], [0.9, 3/12]]
r4 = [[0, 3/12], [0.2, 3/12], [0.4, 3/12], [0.6, 3/12], [0.8, 3/12], [1.0, 3/12]]

v5 = [[0, 4/12], [0.2, 4/12], [0.4, 4/12], [0.6, 4/12], [0.8, 4/12], [1.0, 4/12]]
r5 = [[0.1, 4/12], [0.3, 4/12], [0.5, 4/12], [0.7, 4/12], [0.9, 4/12]]

v = v1 + v2 + v3 + v4 + v5
r = r1 + r2 + r3 + r4 + r5
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

V = FunctionSpace(mesh, 'CG', 1)
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


rho2 = Function(V)  #Structural materials
rho3 = Function(V)  #Responsive materials
rho = Function(V)
rho2.dat.data[:] = rho2_array
rho3.dat.data[:] = rho3_array

rho = rho3 - rho2
rho = interpolate(rho, V)

File("initial/rho.pvd").write(rho)
File("initial/rho2.pvd").write(rho2)
# File("initial/rho2.pvd").write(rho2)
File("initial/rho3.pvd").write(rho3)
