'''
The purpose of this script is to plot the mean of the particle
distribution as a function of time.

Data needs to be of the form (with no header line):
avg_x, avg_x2, avg_x2-avg_x*avg_x
'''
from sys import argv
import matplotlib.pyplot as plt
import numpy as np

scriptname, filename = argv #unpacking the argument variables

with open(filename,'r') as data:
	mean_x = []
	for line in data:
		ls = line.split()
		mean_x += [float(ls[0])]
	ts = len(mean_x)

times = np.array(range(0,ts))
plt.plot(times,mean_x)
plt.show()