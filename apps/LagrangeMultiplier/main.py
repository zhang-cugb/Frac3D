import sys
sys.path.insert(0, '../..')

from dolfin import *
from multiphenics import *
from libs.xdmf2dolfin import *
from libs.dolfin2multiphenics import *
parameters["ghost_mode"] = "shared_facet"


directory = "data"
file = "lowDimFrac3D"

print("Importing mesh to Dolfin.")
msh, boundaries, subdomains = DolfinReader(directory, file)
print("Creating upper subdomain restriction.")
upper = SubdomainRestriction(msh, subdomains, 19)
print("Creating lower subdomain restriction.")
lower = SubdomainRestriction(msh, subdomains, 20)
print("Creating bottom subdomain restriction.")
bottom = SubdomainRestriction(msh, subdomains, 18)
print("Creating upper/lower interface restriction.")
fracture = InterfaceRestriction(msh, subdomains, {19, 20})
print("Creating lower/bottom interface restriction.")
interface = InterfaceRestriction(msh, subdomains, {18, 20})

print("Creating function spaces.")
V = FunctionSpace(msh, "Lagrange", 2)
print("Defining block function spaces.")
W = BlockFunctionSpace([V, V, V, V, V], restrict=[upper, lower, bottom, fracture, interface])
print("Defining trial functions.")
p1p2p3l1l2 = BlockTrialFunction(W)
(p1, p2, p3, l1, l2) = block_split(p1p2p3l1l2)
print("Defining test functions.")
v1v2v3m1m2 = BlockTestFunction(W)
(v1, v2, v3, m1, m2) = block_split(v1v2v3m1m2)

print("Defining subdomains discretized volumes.")
dx = Measure("dx")(subdomain_data=subdomains)
print("Defining boundaries discretized areas.")
ds = Measure("ds")(subdomain_data=boundaries)
dS = Measure("dS")(subdomain_data=boundaries)
dx1, dx2, dx3, dS1, dS2 = dx(19), dx(20), dx(18), dS(29), dS(30)

K1, K2, K3 = 1e-6, 1e-6, 1e-5

print("Defining coefficients matrix blocks.")
A11 = inner(K1*grad(p1), grad(v1))*dx1
A14 = l1("-")*v1("-")*dS1
A22 = inner(K2*grad(p2), grad(v2))*dx2
A24 = -l1("+")*v2("+")*dS1
A25 = l2("-")*v2("-")*dS2
A33 = inner(K3*grad(p3), grad(v3))*dx3
A35 = -l2("+")*v3("+")*dS2
A41 = m1("-")*p1("-")*dS1
A42 = -m1("+")*p2("+")*dS1
A52 = m2("-")*p2("-")*dS2
A53 = -m2("+")*p3("+")*dS2
a = [[A11,	0,		0,		A14,	0],
	 [0,	A22,	0,		A24,	A25],
	 [0,	0,		A33,	0,		A35],
	 [A41,	A42,	0,		0,		0],
	 [0,	A52,	A53,	0,		0]]
print("Assemblying coefficients matrix A.")
A = block_assemble(a)

print("Defining independent vector blocks.")
F1 = v1*dx1
F2 = v2*dx2
F3 = v3*dx3
f = [F1,	F2,		F3,		0,		0]
print("Assemblying independent vector F.")
F = block_assemble(f)

print("Defining Dirichlet boundary conditions.")
bc1 = DirichletBC(W.sub(0), Constant(4.), boundaries, 32)
bc2 = None
bc3 = DirichletBC(W.sub(2), Constant(1.), boundaries, 31)
bc4 = None
bc5 = None
bcs = BlockDirichletBC([bc1, bc2, bc3, bc4, bc5])
print("Assigning Dirichlet boundary conditions to coefficients matrix A.")
bcs.apply(A)
print("Assigning Dirichlet boundary conditions to independent vector F.")
bcs.apply(F)

print("Defining solution vector P.")
P = BlockFunction(W)
print("Solving linear system AP=F.")
block_solve(A, P.block_vector(), F, "mumps")

print("Exporting splitted solution.")
p1_h, p2_h, p3_h, l1_h, l2_h = block_split(P)
p1_h.rename("p1", "p1")
p2_h.rename("p2", "p2")
p3_h.rename("p3", "p3")

print("Writing solution to output file.")
directory = "results"
file = "lowDimFrac3D_solution_splitted"
XDMFFieldWriter(directory, file, [p1_h, p2_h, p3_h])