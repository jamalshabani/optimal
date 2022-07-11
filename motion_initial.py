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
	if (0.02 < s < 0.04):
		return (cos((s - 0.04)*pi/0.02) + 1)/2
	if (s > 0.04):
		return 1

def g_r(s):
	if (s < 0.02):
		return 1
	if (0.02 < s < 0.04):
		return (cos((s - 0.02)*pi/0.02) + 1)/2
	if (s > 0.04):
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
rho2 = Function(V)  # Structural material 1(Blue)
rho3 = Function(V)  # Responsive material 2(Red)

rho2.dat.data[:] = rho2_array
rho3.dat.data[:] = rho3_array

rho2 = Constant(0.5)
rho3 = Constant(0.5)
