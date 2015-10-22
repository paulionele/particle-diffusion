import numpy as np
import matplotlib.pyplot as plt

with open('diffusion_coeffs_v_probablities.txt','r') as f1:
	p = []
	D = []
	for line in f1:
		ls = line.split()
		p += [float(ls[0])]
		D += [float(ls[1])]

plt.plot(p,D,'-ro')
plt.xlabel('Probability')
plt.ylabel('Diffusion Coefficient')
plt.grid('on')
plt.show()