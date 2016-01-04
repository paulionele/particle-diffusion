'''
Plotter script for 1d_diffusion. Replaced by ipython notebook.
This script is otherwise the basis for a future plotter script
if I can get it to work.

Known issues: 
1) Newly generated plots do not display without having to close 
current figure window. Ideally, the window would remain open and
the figure axis would persist, only redrawing the plots and scaling
as necessary.
2) Not vectorized using numpy. The addition of a second plot such
as a Fokker-Plank distribution is easiest if arrays are vectorized.
'''

from os import path
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import time as ti

#Storing the times and distributions as lists, if that's useful...
times = 0
dists = []

savepath = '/home/paul/Documents/thesis/particle-diffusion/data/'
filename = 'FILE_NAME_HERE.txt'
filepath = path.join(savepath,filename)
with open(filepath,'r') as f:
    
    dl = f.readline()
    dls = list(dl.split())
    dists += [[int(x) for x in dls]]

    times += 1
    number_cells = len(dls)
    number_parti = max(dls) #number of particles in the distribution
    number_index = dls.index(number_parti) #returns an index for the location of start

    for dl in f:
        times += 1
        dls = list(dl.split())
        dists += [[int(x) for x in dls]]

cells = [x for x in range(0,number_cells)]

for time in [0,1,3,10,400,1000,3000,4999]:
    plt.plot(cells,dists[time])
    plt.title('time: {}'.format(time))
    plt.show()

'''
What is returned has the same form as below.
All values are either floats or integers.
The dimensions of a dist should be equal to that of cells.
cells = [0,1,2,3,4,..., number_cells]
times = [t1,t2,t3, ... ]
dists = [[1,2,3,4], [2,3,4,5], [3,4,5,6], ... ]
'''
