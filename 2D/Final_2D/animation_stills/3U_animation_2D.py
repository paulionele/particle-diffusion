'''
Script for plotting 2D density distribution at a selected time step.
Compared to the 1D script, here we do not overlay MC and ME computations.

Input data form:
Each 2D array is separated by an empty line. This empty line is tested for and used to
indirectly count the number of 2D arrays in the data file to determine the number of time steps.

This empty line also defines the start of a new 2D array. For each 2D array:
    -each new line is a new row of data (dimension m).
    -the number of elements in each row (dimension n).
    -for each array, a list of lists is created. Each inner list is a row of data and the whole list
    is contains all data for a single array.
    -The data structure is as follows:
    total_array = [...,array_ith,...] = [...,[...,[array_ith_row_m],...],...]
'''

from os import path
import numpy as np
import matplotlib.pyplot as plt

'''Parsing the input file.'''

times = 0

savepath = '/home/paul/Documents/thesis/particle-diffusion/2D/Final_2D/animation_stills/animation_data_3U/'
#filename = 'AFD_t-20k_xC-15_xE-15_yC-15_yE-15_nU-3_pi-0.05_pe-0.2_pie-0.025.txt'
filename = 'AMC_t-12k_N-500k_nU-3_pi-0.05_pe-0.2_pie-0.025.txt'
filepath = path.join(savepath,filename)

with open(filepath,'r') as f:
  #Open file and count empty lines until desired line is reached.
  # After, record all non-empty lines until empty line is reached again.
  array_ith = []
  time_capture = 10500
  t = 0

  for dl in f:
    if (dl != '\n' and t == time_capture):
      #If the data line is not empty and this is the array to capture...
      dls = [float(x) for x in dl.split()]
      array_ith += [dls]
    elif dl == '\n':
      #The line is empty, increment the counter.
      t += 1
    else:
      pass

plt.imshow(array_ith, cmap='jet', interpolation='nearest')
plt.savefig(savepath + 'heterogeneous_3U_2D.pdf')
plt.show()