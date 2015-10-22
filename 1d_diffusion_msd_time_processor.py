'''
The purpose of this script is to determine to diffusion coefficient for a
particular forward step (unbiased) probability. The calculated diffusion
coeffcient is returned and written to a file along with its cooresponding 
probability. A different script is used to plot that relation.

Data needs to be of the form (with no header line):
avg_x, avg_x2, avg_x2-avg_x*avg_x
The last difference is referred to as msd_values.
'''
from sys import argv
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np

scriptname, filename = argv #unpacking the argument variables

with open(filename,'r') as f1:
	#mean_x = []
	#mean_x2 = []
	msd_values = []
	for line in f1:
		ls = line.split()
		#mean_x += [float(ls[0])]
		#mean_x2 += [float(ls[1])]
		msd_values += [float(ls[2])]
	ts = len(msd_values)

times = np.array(range(0,ts))
regression_information = stats.linregress(times,msd_values)
(slope, intercept, r_value, p_value, std_error) = regression_information
diffusion_coefficient  = slope/2

print('The slope: {}'.format(slope))
print('The diffusion coeff: {}'.format(diffusion_coefficient))

#Writing probability and diffusion coefficient data to file...
with open('diffusion_coeffs_v_probablities.txt','a') as f2:
	
	#Searching for probability in file name. Probability of form 'xyz'.
	p = ''
	a = []
	for i in filename:
		if i.isdigit():
			a += [i]
	for j in a:
		p += j
	probability = int(p[2:4])/100
	f2.write('{} {}\n'.format(probability, diffusion_coefficient))