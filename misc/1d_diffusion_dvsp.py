import numpy as np
import matplotlib.pyplot as plt

def function(a,t):
	p = np.linspace(0.0,0.5)
	D = a*a*p/t
	return p, D

with open('diffusion_coeffs_v_probablities.txt','r') as f1:
	p = []
	D = []
	for line in f1:
		ls = line.split()
		p += [float(ls[0])]
		D += [float(ls[1])]

plt.plot(p,D,'-ro')
plt.plot(*function(1,1))
plt.xlabel('Probability')
plt.ylabel('Diffusion Coefficient')
plt.grid('on')
plt.show()