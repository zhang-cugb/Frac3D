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
print("Creating upper/lower fracture restriction.")
fracture = InterfaceRestriction(msh, subdomains, {19, 20})
print("Creating lower/bottom interface restriction.")
interface = InterfaceRestriction(msh, subdomains, {18, 20})

print("Defining function spaces.")
V = FunctionSpace(msh, "CG", 2)
Gamma = FunctionSpace(msh, "DG", 0)
print("Creating mixed function space.")
W = BlockFunctionSpace([V, V, V, V, V], restrict=[upper, fracture, lower, bottom, interface])
print("Defining trial functions.")
p31p2p32p33l = BlockTrialFunction(W)
(p31, p2, p32, p33, l) = block_split(p31p2p32p33l)
print("Defining test functions.")
v31v2v32v33m = BlockTestFunction(W)
(v31, v2, v32, v33, m) = block_split(v31v2v32v33m)

print("Calculating subdomains discretized volumes.")
dx = Measure("dx")(subdomain_data=subdomains)
dx31, dx2, dx32, dx33 = dx(19), dx(29), dx(20), dx(18)
print("Calculating boundaries and interfaces discretized areas.")
ds = Measure("ds")(subdomain_data=boundaries)
dS = Measure("dS")(subdomain_data=boundaries)
dS12, dS23 = dS(29), dS(30)

print("Assigning subdomains permeability.")
K_values = {}
K_values[19] = Constant(1e-6)
K_values[29] = Constant(1e-3)
K_values[20] = Constant(1e-6)
K_values[18] = Constant(1e-5)
K = Function(Gamma)
for cell in range(len(subdomains.array())):
	subdomain = subdomains.array()[cell]
	K.vector()[cell] = K_values[subdomain]
kappa = Constant(20)

print("Defining coefficients matrix blocks.")
A11 = K*inner(grad(p31), grad(v31))*dx31
A21 = -kappa*p31("-")*v2("-")*dS12
A22 = K*inner(grad(p2), grad(v2))*dx2 + kappa*p2("+")*v2("+")*dS12 + kappa*p2("-")*v2("-")*dS12
A23 = -kappa*p32("+")*v2("+")*dS12
A33 = K*inner(grad(p32), grad(v32))*dx32
A35 = l("-")*v32("-")*dS23
A44 = K*inner(grad(p33), grad(v33))*dx33
A45 = -l("+")*v33("+")*dS23
A53 = m("-")*p32("-")*dS23
A54 = -m("+")*p33("+")*dS23
a = [[A11,	0,		0,		0,		0],
	 [A21,	A22,	A23,	0,		0],
	 [0,	0,		A33,	0,		A35],
	 [0,	0,		0,		A44,	A45],
	 [0,	0,		A53,	A54,	0]]
print("Assemblying coefficients matrix A.")
A = block_assemble(a)

print("Defining independent vector blocks.")
F1 = Constant(0.)*v31*dx31
f = [F1,	0,		0,		0,		0]
print("Assemblying independent vector F.")
F = block_assemble(f)

print("Defining Dirichlet boundary conditions.")
bc1 = DirichletBC(W.sub(0), Constant(4.), boundaries, 32)
bc2 = None
bc3 = None
bc4 = DirichletBC(W.sub(3), Constant(1.), boundaries, 31)
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
p31_h, p2_h, p32_h, p33_h, l_h = block_split(P)
p31_h.rename("p31", "p31")
p2_h.rename("p2", "p2")
p32_h.rename("p32", "p32")
p33_h.rename("p33", "p33")

print("Writing solution to output file.")
directory = "results"
file = "lowDimFrac3D_solution"
XDMFFieldWriter(directory, file, [p31_h, p2_h, p32_h, p33_h])