import matplotlib.pyplot as plt
#from time import sleep
#import matplotlib.mlab as mlab
import matplotlib.animation as animation
import numpy as np 

'''
Need to work on file selection...
The following code produces a list of times, a list of lists of particle
distributions, and an array of values cooresponding to cell numbers for plotting.
'''

times = []
dists = []
filename = 'particles-10000_xcells-80_time-100_startpos-40.txt'
with open(filename,'r') as f:

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

# def setup_backend(backend='TkAgg'):
#     import sys
#     del sys.modules['matplotlib.backends']
#     del sys.modules['matplotlib.pyplot']
#     import matplotlib as mpl
#     mpl.use(backend)  # do this before importing pyplot
#     import matplotlib.pyplot as plt
#     return plt

def animate():
    mu, sigma = 100, 15
    N = 4
    x = mu + sigma * np.random.randn(N)
    rects = plt.bar(range(N), x, align='center')
    for i in range(50):
        x = mu + sigma * np.random.randn(N)
        for rect, h in zip(rects, x):
            rect.set_height(h)
        fig.canvas.draw()

#plt = setup_backend()
fig = plt.figure()
win = fig.canvas.manager.window
win.after(10, animate)
plt.show()


#x and y should have the same length.
# for time in range(0,len(times)-98):
# 	distribution = dists[time]
# 	plt.bar(cells,distribution, facecolor='green',alpha=0.5, align='center')
# 	plt.xlabel('Cells')
# 	plt.ylabel('Particles')
# 	plt.title('Time = {0}'.format(str(time)))
# 	plt.show()
# 	plt.clf()

	# add a 'best fit' line
	#y = mlab.normpdf( bins, distribution_mean, distribution_width)
	#l = plt.plot(bins, y, 'r--', linewidth=1)