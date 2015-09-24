import random
#import sys

#Defined constants
N = 1000       #number of particles
x_cells = 21    #number of cells in x (just a 1D array for now)
end_time = 100  	#number of time steps

steps = [t for t in range(0,end_time)]
cells = [0*n for n in range(0,x_cells)] #1D array of cells, elements intilized to zero
#if x_cells % 2 != 0
#	start_pos1 = int(x_cells/2) + 1
#	cells[start_pos1] = N
#else:
#	start_pos1 = int(x_cells/2)
#	cells[start_pos1]  = N
start_pos1 = 11
cells[start_pos1] = N

print('Dimensions of array: {m}x{n}'.format(m = x_cells, n = 1))
print('The number of time steps: {t}'.format(t = len(steps)))
print('The number of particles: {n}'.format(n = N))
'''
As for the structure to hold the data- can be flexible with this. I've decided 
to put everything into a list. Each element of the list will be a 2-tuple. The
first object will hold the time step, the second object a list of values of the
current distribution of the particles (basically the population of each cell at
that particular time step).

data = [(t1,[0,0,0,10,0,...]),(t2,[0,0,2,5,3,0,...]),...]
'''
data = [(-1,cells)] #need to shift indicies over.

for time in steps:
	'''
	For each time step, compute over all cells.
	Create an empty array for the particles to be added into during randomized
	movement. Once all the particles have been moved, the old distribution is 
	replaced with the new one.
	'''
	cells_current = [0*n for n in range(0,x_cells)]
	for i in range(0,x_cells):
		'''
		The index for the cells is 'i'. Note indexing starts at zero.
		The index ends at (x_cells - 1), (nature of the range function).
		Basically- for each cell- we'll loop over all the particles and randomly
		decide if a particle moves. If it does move, then randomly decide if to 
		the left or to the right. Particles at the first and last cell are 
		subject to boundary conditions, again a random choice decides if the 
		particle stays in place or moves away from the boundary cell.

		Also considering empty cells. Cannot subtract from empty cell. The next
		if suite tests for this. The next iteration of the loop begins if the
		i'th cell is found empty.

		Do not want particles to be handled twice in one time step. I've written
		the program in such a way for each particle movement iteration, a 
		current particle distribution is updated. Basically this preserves
		the 'old' distribution while particles are moved and then replaces
		that distribution.
		'''
		if not cells[i]:
			#go back to the start of the loop if cell is empty
			continue

		for particle in range(0,cells[i]):
			'''
			Processes all the particles in the i'th cell.
			Particles either moved left, right, or stay in current cell. 
			'''
			if random.random() < 0.5:
				#particle will move, next draw decides direction.
				lor = random.random()
				if lor < 0.5 and i != 0 and i != (x_cells - 1): #take this -1 away?):
					#particle will move left, +1 left cell
					cells_current[i-1] += 1
				elif lor > 0.5 and i != 0 and i != (x_cells - 1): #take this -1 away?
					#particle will move right, +1 right cell
					cells_current[i+1] += 1
				
				#the following are the boundary cases
				elif i == 0:
					if lor < 0.5:
						#particle moves to right
						cells_current[i+1] += 1
					else:
						#particle stays in place
						cells_current[i] += 1
				elif i == (x_cells - 1):
					if lor < 0.5:
						#particle moves to left
						cells_current[i-1] += 1
					else:
						#particle stays in place
						cells_current[i] += 1
				else:
					print('error')
			else:
				#particle will not move
				cells_current[i] += 1
	cells = cells_current
	data += [(time,cells)] #10% increase in speed

#writing data to file...
filename = 'particles-{0}_xcells-{1}_time-{2}_startpos-{3}.txt'.format(
			N,x_cells,end_time,start_pos1)	
with open(filename,'w') as f:
	for time_step_data in data:
		f.write(str(time_step_data)+'\n')
print('done')
