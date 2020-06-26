from dolfin import Mesh, MeshValueCollection, XDMFFile
from dolfin.cpp.mesh import MeshFunctionSizet, SubsetIterator


def DolfinReader(directory, file):
	msh = Mesh()
	with XDMFFile("{}/{}.xdmf".format(directory, file)) as infile:
		infile.read(msh)
		dim = msh.topology().dim()
	mvc = MeshValueCollection("size_t", msh, dim=dim-1)
	with XDMFFile("{}/{}_facets.xdmf".format(directory, file)) as infile:
		infile.read(mvc, "boundaries")
	boundaries = MeshFunctionSizet(msh, mvc)
	mvc = MeshValueCollection("size_t", msh, dim=dim)
	with XDMFFile("{}/{}_physical_region.xdmf".format(directory, file)) as infile:
		infile.read(mvc, "subdomains")
	subdomains = MeshFunctionSizet(msh, mvc)
	return msh, boundaries, subdomains


def XDMFFieldWriter(directory, file, solutionList):
	with XDMFFile("{}/{}.xdmf".format(directory, file)) as infile:
		infile.parameters["rewrite_function_mesh"] = False
		infile.parameters["functions_share_mesh"] = True
		for solution in solutionList:
			infile.write(solution, 0.0)
		infile.close()
	return


def JoinSolutions(solutionList, subdomains, subdomainsList, V):
	dofmap = V.dofmap()
	mesh = V.mesh()
	joinedSolution = solutionList[1]
	for i,ID in enumerate(subdomainsList):
		cells = SubsetIterator(subdomains, ID)
		for cell in cells:
			dofs = dofmap.cell_dofs(cell.index())
			for dof in dofs:
				joinedSolution.vector()[dof] = solutionList[i].vector()[dof]
	return joinedSolution