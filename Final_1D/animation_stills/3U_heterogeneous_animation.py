'''
This is a plotting script specific to a set of data, in particular the heterogeneous 1D
simulation data with square regions of 15 lattice sites and a total of 11 unit cells.

Plots produced are layered onto a single figure. MC, ME, and FP plots on single figure.
'''

import os
import numpy as np
import matplotlib.pyplot as plt

def FP(time, D, cells, number_index):
  '''
  This function generates a Fokker-Planck distribution for selected time step.
  Output is an array of length equal to the total array length. Plot against 'cells'.
  ---Need to calculate an effective D in the long time using MSD plot.
  '''
  P = 1.0/(np.sqrt(4.0*np.pi*D*float(time))) * np.exp(-((cells-number_index)*(cells-number_index))/(4.0*D*float(time)))
  return number_parti*P

def extractor(savepath, filename):
  '''
  Function called to open file, extract the density distribution data from
  the desired time step and return the time step and distribution data for 
  plotting. Modifications are applied to data before returning.
  '''
  #Need filename to apply appropriate modifications to distribution.
  filepath = os.path.join(savepath, filename)
  time_capture = 19000
  distribution = []

  with open(filepath,'r') as f:

    for i, line in enumerate(f):
      #Pick out specific lines in the selected file.
      if i == 0:
        #Extraction of basic information.
        ls = [float(x) for x in line.split()]
        length = len(ls) #number of lattice points in the distribution
        index = ls.index(max(ls)) #returns an index for the location of start (CAUTION)

      elif i == (time_capture + 1):
        #Extraction of density distribution data.
        ls = [float(x) for x in line.split()]
        dd = np.array(ls)
      else:
        pass

  return length, dd, index

N = 800000
savepath = '/home/paul/Documents/thesis/particle-diffusion/Final_1D/animation_stills/heterogeneous_data_3U/'
files = sorted(os.listdir(savepath))



for filename in files:
  if '.txt' in filename:
    #Get only text files.
    (length, dd, index) = extractor(savepath, filename)
    cells = np.array([x for x in range(0,length)])
    
    #Plot some region-defining lines first so the don't overlay other plots.
    for i in range(0,length):
      if i%15 == 0:
        plt.plot((i, i), (0, 1/25), color = 'silver', linestyle = '--')

    if 'MC' in filename:
      #For the MC plot, data normalization and custom labels.
      plt.plot(cells, (1/N)*dd, label = 'Monte Carlo')
      plt.xlabel('Lattice Site', fontsize = 16)
      plt.ylabel('Density', fontsize = 16)
      plt.xlim(0,length-1)
      plt.ylim(0,1/33)
      plt.legend(loc='upper right')

    if 'FD' in filename:
      plt.plot(cells, dd, label = 'Master Equation')
      plt.xlabel('Lattice Site', fontsize = 16)
      plt.ylabel('Density', fontsize = 16)
      plt.xlim(0,length-1)
      plt.ylim(0,1/33)
      plt.legend(loc='upper right')

  else:
    print('Non *.txt hit.')

plt.savefig(savepath + '3U_heterogeneous_plots_1D.pdf')
plt.show()


#plot fokker-plank distribution DOES NOT WORK PROPERLY ANYMORE :(
#     if time != 0:
#         plt.plot(cells, FP(time, 0.1, cells, number_index))    