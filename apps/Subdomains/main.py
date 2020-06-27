import sys
sys.path.insert(0, '../..')

from dolfin import *
import numpy
from libs.xdmf2dolfin import *
parameters["ghost_mode"] = "shared_facet"


directory = "data"
file = "lowDimFrac3D"

print("Importing mesh to Dolfin.")
msh, boundaries, subdomains = DolfinReader(directory, file)

print("Creating function spaces.")
V = FunctionSpace(msh, "CG", 2)
K = FunctionSpace(msh, "DG", 0)
print("Defining trial functions.")
p = TrialFunction(V)
print("Defining test functions.")
v = TestFunction(V)

print("Defining subdomains discretized volumes.")
dx = Measure("dx")(subdomain_data=subdomains)
print("Defining boundaries discretized areas.")
ds = Measure("ds")(subdomain_data=boundaries)

print("Assigning subdomains permeability.")
k_values = {}
k_values[19] = 1e-6
k_values[20] = 1e-6
k_values[18] = 1e-5
k = Function(K)
for cell in range(len(subdomains.array())):
	subdomain = subdomains.array()[cell]
	k.vector()[cell] = k_values[subdomain]

print("Defining bilinear form.")
a = k*inner(grad(p), grad(v))*dx

print("Defining linear form.")
L = Constant(0.)*v*dx

print("Assigning Dirichlet boundary conditions.")
bc1 = DirichletBC(V, Constant(4.), boundaries, 32)
bc2 = DirichletBC(V, Constant(1.), boundaries, 31)
bcs = [bc1, bc2]

print("Creating variational problem a(p,v)=L(p).")
P = Function(V)
problem = LinearVariationalProblem(a, L, P, bcs)
solver = LinearVariationalSolver(problem)
solver.parameters["linear_solver"] = "mumps"
solver.solve()
P.rename("p", "p")

print("Writing solution to output file.")
directory = "results"
file = "lowDimFrac3D_solution"
XDMFFieldWriter(directory, file, [P])