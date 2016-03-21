'''
The purpose of this script is to create plots of the effective diffusion
coefficient vs. time, and overlay them all onto a single graph.

INPUT: Multiple files, selected automatically from folder using utility.
OUTPUT: A single graph with multiple superposed plots.

Can probably vectorize some lists to make everything faster.
'''
import os
import matplotlib.pyplot as plt
import numpy as np

def huh(savepath, filename):
	#Call to process a file's contents.
	filepath = os.path.join(savepath,filename)
	with open(filepath,'r') as f1:
		#col1 = [] #typically mean x-position
		#col2 = [] #typically mean (x-position)^2
		col3 = [] #typically MSD-x as calculated from the above quantities
		for line in f1:
			ls = line.split()
			#Analytical data output is currently set at 3 cols.
			#Make sure you are using the correct data. 
			#col1 += [float(ls[0])]
			#col2 += [float(ls[1])]
			col3 += [float(ls[2])]

	times = np.array(range(1, len(col3) + 1)) #create an array of times
	msd_x = np.array(col3) #create an array of msd_x
	deff  = np.true_divide(msd_x, 2*times)

	return times, deff

#Getting a list of file names in current directory (caution!).
directory = 'c_yE/'
savepath = '/home/paul/Documents/thesis/particle-diffusion/2D/2D-data/' + directory
files = sorted(os.listdir(savepath))

#Commands for latex labelling.
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

for filename in files:
	if '.txt' in filename:
		#Get only text files.
		(times, deff) = huh(savepath, filename)
		plot_name = filename.split('_')[1]

		plt.plot(times, deff, label = plot_name)
		plt.xscale('log')
		plt.yscale('log')
		plt.xlabel('Time Step', fontsize = 16)
		plt.ylabel('Effective Diffusion Coefficient', fontsize = 16)
	else:
		pass

plt.legend(loc=2)
plt.savefig(savepath + 'deff_time.pdf')
plt.show()