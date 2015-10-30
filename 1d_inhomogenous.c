/* 
Diffusion simulation on 1D, inhomogenous lattice.

The unit cell is built so that the overall cell lattice is non-symmetric. That's
something I need to fix.. I like symmetry.
*/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

// CPP directives.
#define N 1000  //number of particles
#define xC 20    //x length of cellular space
#define xE 5    //x length of extracellular space
#define nU 8  //x number of unit cells

int main(){
	//x unit cell length (xU), total array length (xL) and centered start (SP).
	int xU = xC + xE;
	int xL = xU*nU;
	int SP = (int)(xL/2.0);
	int modx;

	//Index variables for the for-loops and time-step limit.
	int i;
	int t, tmax = 30;

	//Stepsize and probabilities for intracellular (0) and extracellular diffusion (1).
	//Probablities crossing from intra to extra (pie) or extra to intra (pei). 
	float a = 1.0;
	float pli = 0.3, pri = 0.3, psi = 1.0-pli-pri;
	float ple = 0.2, pre = 0.2, pse = 1.0-ple-pre;
	float pie = 0.1, pei = 0.1; 

	//Random variable and variables for MSD calculations.
	float rnd;
	float sum_x, sum_x2, avg_x, avg_x2;

	//Shifted index (particle pos. + SP) and an array to hold particle positions.
	int sxi;
	float x[N] = {0};

	//Density distribution (rho) and setting the start position.
	int *rho;
	rho = malloc(xL*sizeof(int));
	rho[(int)SP-1] = N;

	//Some file naming and preparation.
	char *f1 = "/home/paul/Documents/thesis/particle-diffusion/data/\
		TEST1.txt";
	char *f2 = "/home/paul/Documents/thesis/particle-diffusion/data/\
		TEST1-stats.txt";
	FILE *dd, *outstats;
	dd = fopen(f1, "w");
	outstats = fopen(f2, "w");
	
	//Writing in the initial ditribution.
	fprintf(dd, "#First integer is the time step.#\n");
	for(i = 0; i < xL; i++){fprintf(dd, "%d ", rho[i]);}
	fprintf(dd, "\n");

	//--------------------------------------------------------------------------
	//Below is where the main work is done.
	for (t = 0; t < tmax; t++){
		printf("%d\n",t);
		sum_x = 0;
		sum_x2 = 0;

		//Loop over every cell in array, set particle count to 0.
		for (i = 0; i < xL; i++){rho[i] = 0;}

		for(i = 0; i < N; i++){
			//Loop over all particles. sxi: shifted position for ith particle.
			//The modx value returned is 0:(xU-1). This can be used to determine
			// if the particle is at a boundary, i.e if modx = 0 or (xU-1).
			sxi = (int)(x[i] + SP);
			modx = sxi % xU;

			//The order of the unit cell can be reversed from here.
			if(modx < xC){
				//Perform operations specific to INTRACELLULAR conditions.

				//Check the ith particle position and determine if at boundary.
				if(modx != 0 && modx != (xU-1)){
					//Particle is not at inner boundary, may move left, right, or stay.
					rnd = (float)rand()/(float)RAND_MAX;
					if(rnd < pli){
						//The ith particle moves a distance 'a' to the left.
						x[i] -= a;
					}
					else if(rnd > pli && rnd < (pli + pri)){
						//The ith particle moves a distance 'a' to the right.
						x[i] += a;
					}
				}

				//Intracellular boundary conditions.
				else if(modx == 0){
					rnd = (float)rand()/(float)RAND_MAX;
					if(rnd < pie){
						//The ith particle moves a distance 'a' to the left and through boundary.
						x[i] -= a;
					}
				}
				else if(modx == (xU-1)){
					rnd = (float)rand()/(float)RAND_MAX;
					if(rnd < pie){
						//The ith particle moves a distance 'a' to the right and through boundary.
						x[i] += a;
					}
				}
			}
			else{
				//Perform operations specific to EXTRACELLULAR conditions.

				//Check the ith particle position and determine if at boundary.
				if(modx != 0 && modx != (xU-1)){
					//Particle is not at boundary, may move left, right, or stay.
					rnd = (float)rand()/(float)RAND_MAX;
					if(rnd < ple){
						//The ith particle moves a distance 'a' to the left.
						x[i] -= a;
					}
					else if(rnd > ple && rnd < (ple + pre)){
						//The ith particle moves a distance 'a' to the right.
						x[i] += a;
					}
				}

				//Extracellular boundary conditions.
				else if(modx == 0){
					rnd = (float)rand()/(float)RAND_MAX;
					if(rnd < pei){
						//The ith particle moves a distance 'a' to the left and through boundary.
						x[i] -= a;
					}
				}
				else if(sxi == (xL-1)){
					rnd = (float)rand()/(float)RAND_MAX;
					if(rnd < pei){
						//The ith particle moves a distance 'a' to the right and through boundary.
						x[i] += a;
					}
				}
			}

			/**/
			sum_x += x[i];
			sum_x2 += x[i]*x[i];

			/*Shift the particle position and increment particle count at that
			position rho[sxi] by 1 unit.*/
			sxi = (int)(x[i] + SP);
			rho[sxi]++; 
		}
		
		/*Write the cell density data to file. First column is time.*/
		fprintf(dd, "%d ",t);
		for(i = 0; i < xL; i++){fprintf(dd, "%d ", rho[i]);}
		fprintf(dd,"\n");
		
		/*Write the mean-square-displacement data to file.*/ 
		avg_x = sum_x/(float)N;
		avg_x2 = sum_x2/(float)N;
		fprintf(outstats, "%f %f %f\n", avg_x, avg_x2, avg_x2-avg_x*avg_x);
	}
}

// get_argv(){
// 	//Get command-line arguements.
// 	//
// }