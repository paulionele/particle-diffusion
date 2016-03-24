'''
This is a plotting script specific to a set of data, in particular the heterogeneous 1D
'pie' data where we plot beta and deff vs. time for varying boundary transition probabilities (pie).

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

def get_beta_deff(savepath, filename):
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

  #For effective diffusion:
  deff  = np.true_divide(msd_x, 2*times)

  #For beta(t) plot.
  log_times = np.log10(times)
  log_msd_x = np.log10(msd_x)
  dt = np.gradient(log_times)
  dd = np.gradient(log_msd_x, dt)

  return times, deff, log_times, dd, msd_x

#Getting a list of file names in current directory (caution!).
savepath = '/home/paul/Documents/thesis/particle-diffusion/Final_1D/analysis_plots/pie_data/'
files = sorted(os.listdir(savepath))

l = 0 #for label incrementation
#labels = ['P_{ie} = 0.005', 'P_{ie} = 0.010', 'P_{ie} = 0.050', 'P_{ie} = 0.100', 'P_{ie} = 0.150']
labels = ['pie = 0.005', 'pie = 0.010', 'pie = 0.050', 'pie = 0.100', 'pie = 0.150']
values = [0.005, 0.010, 0.050, 0.100, 0.150]

f, ax = plt.subplots(2,1)

for filename in files:
  #Files should have been sorted well so label application should be OK.
  if '.txt' in filename:
    #Get only text files.
    (times, deff, log_times, dd, msd_x) = get_beta_deff(savepath, filename)
    
    if 'MC' or 'FD' in filename:
      #For the deff plots.
      if 'MC' in filename:
        #For specific colour and line style.
        plt.sca(ax[1])
        plt.plot(times, deff, linestyle = '--', label = labels[l])
      else:
        plt.sca(ax[1])
        plt.plot(times, deff, label = labels[l])
      
      #Specific marking lines.
      plt.sca(ax[1])
      plt.plot((0, len(times)), (0.05, 0.05), linestyle = '--', color = 'black')
      plt.sca(ax[1])
      plt.plot((0, len(times)), (0.005, 0.005), linestyle = '--', color = 'black')

      plt.xscale('log')
      #plt.yscale('log')
      plt.xlabel('Time', fontsize = 16)
      plt.ylabel('Effective Diffusion Coefficient', fontsize = 16)
      plt.legend(loc = 'center left')

    if 'FD' in filename:
      #For the beta(t) plots.
      plt.sca(ax[0])
      plt.plot((0, np.log10(len(times))), (1.0, 1.0), linestyle = '--', color = 'black')
      plt.sca(ax[0])
      plt.plot(log_times, dd, label = labels[l])
      plt.xlabel('log(Time)', fontsize = 16)
      plt.ylabel('Beta Factor', fontsize = 16)
      plt.legend(loc = 'center left')

    l += 1
  else:
    print('Non *.txt file hit.')

plt.savefig(savepath + 'pie_deff_beta_t_1D.pdf')
plt.show()