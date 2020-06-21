import sys
sys.path.insert(0, '../..')

from dolfin import *
from multiphenics import *
from libs.xdmf2dolfin import *
from libs.dolfin2multiphenics import *


directory = "data"
file = "lowDimFrac3D"

msh, boundaries, subdomains = DolfinReader(directory, file)
upper = SubdomainRestriction(msh, subdomains, 19)
lower = SubdomainRestriction(msh, subdomains, 20)
bottom = SubdomainRestriction(msh, subdomains, 18)
fracture = InterfaceRestriction(msh, subdomains, {19, 20})
interface = InterfaceRestriction(msh, subdomains, {18, 20})

V = FunctionSpace(msh, "Lagrange", 2)
W = BlockFunctionSpace([V, V, V, V, V], restrict=[upper, lower, bottom, fracture, interface])
p1p2p3l1l2 = BlockTrialFunction(W)
(p1, p2, p3, l1, l2) = block_split(p1p2p3l1l2)
v1v2v3m1m2 = BlockTestFunction(W)
(v1, v2, v3, m1, m2) = block_split(v1v2v3m1m2)

dx = Measure("dx")(subdomain_data=subdomains)
ds = Measure("ds")(subdomain_data=boundaries)
dx1, dx2, dx3, ds1, ds2 = dx(19), dx(20), dx(18), ds(29), ds(30)

K1, K2, K3 = 1e-6, 1e-6, 1e-5

A11 = inner(K1*grad(p1), grad(v1))*dx1
A14 = l1("-")*v1("-")*ds1
A22 = inner(K2*grad(p2), grad(v2))*dx2
A24 = -l1("+")*v2("+")*ds1
A25 = l2("-")*v2("-")*ds2
A33 = inner(K3*grad(p3), grad(v3))*dx3
A35 = -l2("+")*v3("+")*ds2
A41 = m1("-")*p1("-")*ds1
A42 = -m1("+")*p2("+")*ds1
A52 = m2("-")*p2("-")*ds2
A53 = -m2("+")*p3("+")*ds2
a = [[A11,	0,		0,		A14,	0],
	 [0,	A22,	0,		A24,	A25],
	 [0,	0,		A33,	0,		A35],
	 [A41,	A42,	0,		0,		0],
	 [0,	A52,	A53,	0,		0]]
A = block_assemble(a)

# F1 = v1*dx1
# F2 = v2*dx2
# F3 = v3*dx3
# f = [F1,	F2,		F3,		0,		0]
# F = block_assemble(f)

# bc1 = DirichletBC(W.sub(0), Constant(4.), boundaries, 32)
# bc2 = None
# bc3 = DirichletBC(W.sub(2), Constant(1.), boundaries, 31)
# bc4 = None
# bc5 = None
# bcs = BlockDirichletBC([bc1, bc2, bc3, bc4, bc5])
# bcs.apply(A)
# bcs.apply(F)

# P = BlockFunction(W)
# block_solve(A, P.block_vector(), F)