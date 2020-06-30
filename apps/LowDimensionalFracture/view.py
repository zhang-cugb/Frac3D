import sys
sys.path.insert(0, '../..')

from libs.tsv2plt import *
from benchmark import compare_my_method


directory = "results"
filename = "PlotOverLine"

data = TSVReader(directory, filename)
arc_length = np.array(data['arc_length'])
p31 = np.array(data['p31'])
p32 = np.array(data['p32'])
p33 = np.array(data['p33'])
pressure = []
for i,length in enumerate(arc_length):
	if p31[i] == p32[i] and p31[i] != 0:
		end31 = i
	if p32[i] == p33[i] and p32[i] !=0:
		end32 = i
for i,length in enumerate(arc_length):
	if i <= end31:
		pressure.append(p31[i])
	elif i > end31 and i <= end32:
		pressure.append(p32[i])
	elif i > end32:
		pressure.append(p33[i])
pressure = np.array(pressure)
case = 1
plot_id = "a"
ref_index = 0
line_or_fracture_id = 0
fig = compare_my_method.evaluate(arc_length, pressure, case, plot_id, ref_index, line_or_fracture_id)
fig.show()