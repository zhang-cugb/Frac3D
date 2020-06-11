import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
from libs.msh2xdmf import *

files = CheckFiles()
ClearOutputFolder()
for f in files:
	msh, msh_facets, msh_physical_region = AnsysReader(f)
	msh.prune()
	XdmfWriter(msh, "", f)
	XdmfWriter(msh_facets, "_facets", f)
	XdmfWriter(msh_physical_region, "_physical_region", f)