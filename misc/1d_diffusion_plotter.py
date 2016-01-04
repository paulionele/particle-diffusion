'''
Plotter script for 1d_diffusion.
'''

from os import path
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import time as ti
def fp(D, time, number_cells, number_parts, number_index):
    '''
    Try and avoid usrdef functions, they are slow as shit...
    or maybe it's the type conversions.
    D is the diffusion coefficient.
    number_parts is the number of particles.
    number_index is the index of the starting cell (ie. the dist mean)
    '''
    x = np.arange(0,number_cells,0.001)
    y = []

    for t in range(1,time+1):
        yc = 1.0/np.sqrt(4.0*np.pi*D*float(t)) * \
            np.exp(-(x-number_index)**2 / (4.0*D*float(t)))
        yc = list(yc)
        y += [yc]

    x = list(x)
    return x, y

#Storing the times and distributions as lists, if that's useful...
times = []
dists = []

savepath = '/home/paul/Documents/thesis/particle-diffusion-data/'
filename = 'TEST2.txt'
fpath = path.join(savepath,filename)
with open(fpath,'r') as f:

    next(f) #move past the preamble to the 'data line' dl
    dl = f.readline() #read the first line after preamble, get length and other info
    dls = list(dl.split())
    number_cells = len(dls) - 1 #don't consider the times
    number_parts = max(dls) #number of particles in the distribution
    number_index = dls.index(number_parts) #returns an index for the location of start
    times += [int(dls[0])] #this is just a list of times
    dists += [[int(x) for x in dls[1:]]] #dists is a list of lists

    for time_step_data in f:
        dls = list(time_step_data.split())
        times += [int(dls[0])]
        dists += [[int(x) for x in dls[1:]]]

cells = [x for x in range(0,number_cells)]

'''
What is returned from the above has the same form as below.
All values are either floats or integers.
The dimensions of a dist should be equal to that of cells.
cells = [0,1,2,3,4,..., number_cells]
times = [t1,t2,t3, ... ]
dists = [[1,2,3,4], [2,3,4,5], [3,4,5,6], ... ]

Below is code for plotting and other analyses.
'''
###
'''
The following creates a list of MSD values, one MSD value for each time step.
Then the MSD values are graphed against time.

The MSD value is calculated as follows: 
sum(xi^2*F(xi)) where F(xi)=(number_particles)/(total_particles)
'''

msd_values = []
for time in times:
    #for each time step, calculate the MSD
    distribution = dists[time]
    s = []
    i = 1
    for x in distribution:
        s += [(i**2)*(x/int(number_parts))]
        i += 1
    msd_values += [sum(s)]

#Getting some information ...
regression_information = stats.linregress(times,msd_values)
(slope, intercept, r_value, p_value, std_error) = regression_information

print('The slope: {}'.format(slope))
print('The intercept: {}'.format(intercept))
print('The r-squared value: {}'.format(r_value**2))
print('The diffusion coeff: {}'.format(slope/2))

plt.plot(times,msd_values)
plt.show()

(x,y) = fp(0.25,max(times),number_cells,1,number_index)
print(len(x))
print(len(y[1]))
plt.plot(x,y[2])
plt.show()
plt.plot(x,y[400])
plt.show()
