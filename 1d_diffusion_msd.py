'''
The purpose of this script is to determine to diffusion coefficient for a
particular forward step (unbiased) probability. The calculated diffusion
coeffcient is returned and written to a file along with its cooresponding 
probability. A different script is used to plot that relation.

This script needs to *-stats.txt file generated for a particular probablity.
The name of the script of interest must be set at the filename variable below.
CL arguments of the form: python3 1d_diffusion_msd_time_processor.py

Data needs to be of the form (with no header line):
avg_x, avg_x2, avg_x2-avg_x*avg_x
The last difference is referred to as msd_values.
'''
from os import path
#from sys import argv
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np

def regression(times, msd_values):
	''' 
	Function called for calculating linear regression.
	Returns a set of parameters used for to construct a linear
	fit plot to the computed data. Possible to 'window'.
	'''
	(slope, intercept, r_value, p_value, std_error) = stats.linregress(times, msd_values)
	diffusion_coefficient  = slope/2 #redundent
	return slope, intercept, diffusion_coefficient

#scriptname, filename = argv #unpacking the argument variables
filename = ''
savepath = '/home/paul/Documents/thesis/particle-diffusion/2D/2D-data/'
filepath = path.join(savepath,filename)

with open(filepath,'r') as f1:
	msd_x = []
	msd_y = []
	msd_values = []
	for line in f1:
		ls = line.split()
		msd_x += [float(ls[0])]
		msd_y += [float(ls[1])]
		msd_values += [float(ls[0])]
	ts = len(msd_values)

times = np.array(range(1,ts+1))

#Regression time window. ts: time_start, th: time_halt
ts1 = 0
th1 = 5

#Plotting time window.
ts2 = 0
th2 = ts

#Fitting a linear curve through the data range specified in the regression window.
(slope, B, D) = regression(times[ts1:th1], msd_values[ts1:th1])
print('The effective diffusion D on [{},{}]: {}'.format(ts1, th1, D))

#Plotting the data.
plt.plot(times[ts2:th2], msd_x[ts2:th2], color = 'blue', label = 'MSD_x')
plt.plot(times[ts2:th2], msd_y[ts2:th2], color = 'green', label = 'MSD_y')
#plt.plot(times[ts2:th2], msd_values[ts2:th2], color = 'red')

#plt.plot(times[ts2:th2], 4*0.2*times[ts2:th2], '--r') #extracellular diffusion
#plt.plot(times[ts2:th2], 4*0.05*times[ts2:th2], '--g') #crossing intra to extra

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Time')
plt.ylabel('MSD')
plt.legend()
plt.show()

#Plot of MSD/t vs. t
#plt.plot(times[ts2+1:th2],np.array(msd_values[ts2+1:th2])/times[ts2+1:th2])
#plt.show()

#Can calculate a diffusion coefficient for every time step. Maybe plot that?

#Writing probability and diffusion coefficient data to file...
with open(savepath + 'TEST_diffusion_coeffs_v_probablities.txt','a') as f2:
	
	#Searching for probability in file name. Probability of form 'xyz'.
	#Comment out the code for inhomogenous case.
	# p = ''
	# a = []
	# for i in filename:
	# 	if i.isdigit():
	# 		a += [i]
	# for j in a:
	# 	p += j
	# probability = int(p[2:4])/100
	# f2.write('{} {}\n'.format(probability, diffusion_coefficient))
	f2.write('{}\n'.format(D))