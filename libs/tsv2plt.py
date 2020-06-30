import csv
import matplotlib.pyplot as plt
import numpy as np


def TSVReader(directory, filename):
	csvdata = []
	with open("{}/{}.tsv".format(directory, filename)) as tsvfile:
		reader = csv.reader(tsvfile, delimiter="\t")
		keys = next(reader)
		numberOfArrays = len(keys)
		for row in reader:
			pushedRow = []
			for elem in row:
				pushedRow.append(float(elem))
			csvdata.append(pushedRow)
	csvdata = np.array(csvdata)
	csvdata = csvdata.transpose()
	data = {}
	for i,key in enumerate(keys):
		data[key] = csvdata[i]
	return data