'''
This script is functioning as expected.
Creates a list of values, elements of either 1 or 0 depending on specifications.
'''

j = 1 #flip the j,k bits to start with cellular or extracellular.
k = 0

lic = 4 #number of cells in intracellular unit cell.
lec = 2 #number of cells in extracellular unit cell.
lsc = 7 #total number of super cells, odd number for symmetry.
ct = [] #generated array

for i in range(0, lsc):
	if j!=0:
		for cell in range(0,lic):
			ct.append(1)
		j = 0
		k = 1
	elif k!=0:
		for cell in range(0,lec):
			ct.append(0)
		j = 1
		k = 0
print(ct)
