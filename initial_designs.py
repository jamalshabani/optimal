# Define function that return different initial designs
from firedrake import *
import numpy as np

def create_initial_design_1_to_3(mesh, V):
	# Initial Design
	mesh_coordinates = mesh.coordinates.dat.data[:]

	# Include initial designs files. 3 different initial designs
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
	rho = interpolate(rho, V)

	return rho


def create_initial_design_1_to_1(mesh, V):
	# Initial Design
	mesh_coordinates = mesh.coordinates.dat.data[:]

	# Include initial designs files. 3 different initial designs
	# Define the points
	r1 = [[0.1, 0], [0.3, 0], [0.5, 0], [0.7, 0], [0.9, 0]]

	r2 = [[0, 0.1], [0.2, 0.1], [0.4, 0.1], [0.6, 0.1], [0.8, 0.1], [1.0, 0.1]]

	r3 = [[0.1, 0.2], [0.3, 0.2], [0.5, 0.2], [0.7, 0.2], [0.9, 0.2]]

	r4 = [[0, 0.3], [0.2, 0.3], [0.4, 0.3], [0.6, 0.3], [0.8, 0.3], [1.0, 0.3]]

	r5 = [[0.1, 0.4], [0.3, 0.4], [0.5, 0.4], [0.7, 0.4], [0.9, 0.4]]

	r6 = [[0, 0.5], [0.2, 0.5], [0.4, 0.5], [0.6, 0.5], [0.8, 0.5], [1.0, 0.5]]

	r7 = [[0.1, 0.6], [0.3, 0.6], [0.5, 0.6], [0.7, 0.6], [0.9, 0.6]]

	r8 = [[0, 0.7], [0.2, 0.7], [0.4, 0.7], [0.6, 0.7], [0.8, 0.7], [1.0, 0.7]]

	r9 = [[0.1, 0.8], [0.3, 0.8], [0.5, 0.8], [0.7, 0.8], [0.9, 0.8]]

	r10 = [[0, 0.9], [0.2, 0.9], [0.4, 0.9], [0.6, 0.9], [0.8, 0.9], [1.0, 0.9]]

	r11 = [[0.1, 1.0], [0.3, 1.0], [0.5, 1.0], [0.7, 1.0], [0.9, 1.0]]

	r = r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 + r9 + r10 + r11

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
	# rho = Constant(0.4)
	rho = interpolate(rho, V)

	return rho


def create_initial_design_1_to_6(mesh, V):
	# Initial Design
	mesh_coordinates = mesh.coordinates.dat.data[:]

	# Include initial designs files. 3 different initial designs
	# Define the points
	r1 = [[0, 0], [0.2, 0], [0.4, 0], [0.6, 0], [0.8, 0], [1.0, 0]]

	r2 = [[0.1, 1/24], [0.3, 1/24], [0.5, 1/24], [0.7, 1/24], [0.9, 1/24]]

	r3 = [[0, 2/24], [0.2, 2/24], [0.4, 2/24], [0.6, 2/24], [0.8, 2/24], [1.0, 2/24]]

	r4 = [[0.1, 3/24], [0.3, 3/24], [0.5, 3/24], [0.7, 3/24], [0.9, 3/24]]

	r5 = [[0, 4/24], [0.2, 4/24], [0.4, 4/24], [0.6, 4/24], [0.8, 4/24], [1.0, 4/24]]

	r = r1 + r2 + r3 + r4 + r5

	def dist(x, y, a, b):
		return sqrt((x - a)**2 + (y - b)**2)

	def g(s):
		if (s < 0.03):
			return 1
		if (0.03 < s < 0.04):
			return (cos((s - 0.03)*pi/0.01) + 1)/2
		if (s > 0.04):
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
	rho = interpolate(rho, V)

	return rho
