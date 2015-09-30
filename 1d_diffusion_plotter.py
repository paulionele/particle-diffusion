from os import path
import matplotlib.pyplot as plt
#from time import sleep
#import matplotlib.mlab as mlab
#import matplotlib.animation as animation
import numpy as np 

'''
Need to work on file selection...
The following code produces a list of times, a list of lists of particle
distributions, and an array of values cooresponding to cell numbers for plotting.
'''

times = []
dists = []

savepath = '/home/paul/Documents/thesis/particle-diffusion-data/'
filename = 'particles-30000_xcells-99_time-500_startpos-50.txt'
fpath = path.join(savepath,filename)
with open(fpath,'r') as f:

	next(f) #move past the preamble to the 'data line' dl
	dl = f.readline() #read the first line after preamble, get length info
	dls = list(dl.split())
	number_cells = len(dls) - 1
	times += [int(dls[0])] #this is just a list of times
	dists += [[int(x) for x in dls[1:]]] #dists is a list of lists

	for time_step_data in f:
		dls = list(time_step_data.split())
		times += [int(dls[0])]
		dists += [[int(x) for x in dls[1:]]]

cells = [x for x in range(0,number_cells)]

'''
cells = [0,1,2,3,4,..., number_cells]
times = [t1,t2,t3, ... ]
dists = [[1,2,3,4], [2,3,4,5], [3,4,5,6], ... ]
The code below is for plotting.
'''

lt = [0,1,2,5,10,20,50,70,90,110,150,200,250,300,350,400,450,500]

for time in lt:
	distribution = dists[time]
	plt.bar(cells,distribution, facecolor='blue', alpha=0.5, align='center')
	plt.xlim(0,number_cells)
	plt.ylim(0,max(dists[0]))
	plt.xlabel('Cells')
	plt.ylabel('Particles')
	plt.title('Distribution of Particles at time-step = {0}'.format(str(time)))
	plt.savefig(savepath + filename[:-4] + '_TIME_{}.png'.format(str(time)))
	plt.close()

	
	# add a 'best fit' line
	# y = mlab.normpdf( bins, distribution_mean, distribution_width)
	# l = plt.plot(bins, y, 'r--', linewidth=1)