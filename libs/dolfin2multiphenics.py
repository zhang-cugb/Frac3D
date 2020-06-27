from dolfin import *
from multiphenics import *


def SubdomainRestriction(mesh, subdomains, subdomain_id):
	D = mesh.topology().dim()
	restriction = MeshRestriction(mesh, None)
	for d in range(D + 1):
		mesh_function_d = MeshFunction("bool", mesh, d)
		mesh_function_d.set_all(False)
		restriction.append(mesh_function_d)
	for c in cells(mesh):
		if subdomains[c] == subdomain_id:
			restriction[D][c] = True
			for d in range(D):
				for e in entities(c, d):
					restriction[d][e] = True
	return restriction


def InterfaceRestriction(mesh, subdomains, subdomain_ids):
	assert isinstance(subdomain_ids, set)
	assert len(subdomain_ids) == 2
	D = mesh.topology().dim()
	restriction = MeshRestriction(mesh, None)
	for d in range(D + 1):
		mesh_function_d = MeshFunction("bool", mesh, d)
		mesh_function_d.set_all(False)
		restriction.append(mesh_function_d)
	for f in facets(mesh):
		subdomains_ids_f = set(subdomains[c] for c in cells(f))
		assert len(subdomains_ids_f) in (1, 2)
		if subdomains_ids_f == subdomain_ids:
			restriction[D - 1][f] = True
			for d in range(D - 1):
				for e in entities(f, d):
					restriction[d][e] = True
	return restriction