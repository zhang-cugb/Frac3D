import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '../..')

from dolfin import *
from multiphenics import *
from numpy import isclose
from libs.xdmf2dolfin import *


directory = "data"
file = "lowDimFrac3D"
msh, boundaries, subdomains = DolfinReader(directory, file)