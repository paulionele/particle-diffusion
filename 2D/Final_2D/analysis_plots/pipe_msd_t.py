'''
This is a plotting script specific to a set of data, in particular the heterogeneous 1D
'pipe' data where we plot MSD vs. time for varying boundary transition probabilities (pie).

Plots produced are layered onto a single figure. Only ME plots on single figure.
'''

import os
import numpy as np
import matplotlib.pyplot as plt


def regression(times, col3):
  ''' 
  Probably not used in this script.
  '''
  (slope, intercept, r_value, p_value, std_error) = stats.linregress(times, col3)
  #diffusion_coefficient  = slope/2 #redundent
  return slope, intercept

def get_msd(savepath, filename):
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

  return times, msd_x

#Getting a list of file names in current directory (caution!).
savepath = '/home/paul/Documents/thesis/particle-diffusion/2D/Final_2D/analysis_plots/pipe_data/'
files = sorted(os.listdir(savepath))

l = 0 #for label incrementation
#labels = ['P_{ie} = 0.005', 'P_{ie} = 0.010', 'P_{ie} = 0.050', 'P_{ie} = 0.100', 'P_{ie} = 0.150']
labels = ['pi/pe = 0.005/0.2', 'pi/pe = 0.010/0.2', 'pi/pe = 0.050/0.2', 'pi/pe = 0.100/0.2', 'pi/pe = 0.150/0.2']

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

#plt.figure(figsize=(6,6))
for filename in files:
  #Files should have been sorted well so label application should be OK.
  if '.txt' in filename:
    #Get only text files.
    (times, msd) = get_msd(savepath, filename)
    
    if 'MC' or 'FD' in filename:
      #For the MSD plots.
      plt.plot(times, msd, label = labels[l])
      plt.xscale('log')
      plt.yscale('log')
      plt.xlabel('Time', fontsize = 16)
      plt.ylabel('MSD', fontsize = 16)
      plt.xlim([0,10**6])
      plt.legend(loc=2)

      l += 1
    else:
      print('Eh.')
  else:
    print('Non *.txt file hit.')

#pei = (pxi/pxe)*pie
#pie = 0.01
#pei = (0.005/0.2)*pie
#p_total = (1/2)*(pie+pei)*(0.005+0.2)*(1/2)

plt.plot(times, 2*0.005*times, color = 'black', linestyle = '--') #lowest intracellular
plt.plot(times, 2*0.1*times, color = 'black', linestyle = '--') #some value???

plt.savefig(savepath + 'pipe_msd_2D.pdf')
plt.show()