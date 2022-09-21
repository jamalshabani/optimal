from firedrake import *# Import gmesh

mesh = Mesh('motion_mesh.msh')

V = TensorFunctionSpace(mesh, 'CG', 1)

Id = Identity(mesh.geometric_dimension()) #Identity tensor
e1 = as_vector((1, 0)) # Direction of responsive material
S =  Id - 2 * outer(e1, e1)

S = interpolate(S ,V)

File('S.pvd').write(S)
print("Done")
