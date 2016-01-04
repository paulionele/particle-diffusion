/* 
Diffusion simulation on 1D, inhomogenous lattice.

Boundary condition tests are handled differently here than the other scripts,
however the program does not work as intended. Maybe return to it later.

The unit cell is built so that the overall cell lattice is non-symmetric. That's
something I need to fix.. I like symmetry.
*/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"

// CPP directives.
#define N 1000  //number of particles
#define xC 20    //x length of cellular space
#define xE 5  //x length of extracellular space
#define nU 3  //x number of unit cells

int main(){
	//x unit cell length (xU), total array length (xL) and centered start (SP).
	int xU = xC + xE;
	int xL = xU*nU;
	int SP = (int)(xL/2.0);
	int modx;

	//Index variables for the for-loops and time-step limit.
	int i;
	int t, tmax = 5000;

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

	//File naming and preparation. Why do we get a segmentation fault below?
	//char *path = "/home/paul/Documents/thesis/particle-diffusion/data/";
	//char *f1 = strcat(path,"TEST1.txt");
	//char *f2 = strcat(path,"TEST1-stats.txt");
	char *f1 = "/home/paul/Documents/thesis/particle-diffusion/data/TEST1.txt";
	char *f2 = "/home/paul/Documents/thesis/particle-diffusion/data/TEST1-stats.txt";
	FILE *outdists, *outstats;
	outdists = fopen(f1, "w");
	outstats = fopen(f2, "w");
	
	//Writing in the initial distribution.
	for(i = 0; i < xL; i++){
		fprintf(outdists, "%d ", rho[i]);
	}
	fprintf(outdists, "\n");

	//--------------------------------------------------------------------------
	//Below is where the main work is done.
	for (t = 1; t < tmax; t++){
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

			//INTRACELLUAR AND NOT AT BOUNDARY
			if(modx != 0 && modx != (xC-1)){
				//Particle is not at inner intracellular boundary, may move left, right, or stay.
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

			//EXTRACELLULAR AND NOT AT BOUNDARY
			else if(modx>= xC && modx!=(xU-1)){
				//Particle is not at extracellular boundary, may move left, right, or stay.
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

			//INTRACELLULAR AND AT LEFT BOUNDARY
			else if(modx == 0){
				rnd = (float)rand()/(float)RAND_MAX;
				if(rnd < pie){
					//The ith particle moves a distance 'a' to the left and through boundary.
					x[i] -= a;
				}
			}

			//INTRACELLULAR AND AT RIGHT BOUNDARY
			else if(modx == (xC-1)){
				rnd = (float)rand()/(float)RAND_MAX;
				if(rnd < pie){
					//The ith particle moves a distance 'a' to the right and through boundary.
					x[i] += a;
				}
			}

			//EXTRACELLULAR AND AT LEFT BOUNDARY
			else if(modx == xC){
				rnd = (float)rand()/(float)RAND_MAX;
				if(rnd < pei){
					//The ith particle moves a distance 'a' to the left and through boundary.
					x[i] -= a;
				}
			}

			//EXTRACELLULAR AND AT RIGHT BOUNDARY
			else if(sxi == (xU-1)){
				rnd = (float)rand()/(float)RAND_MAX;
				if(rnd < pei){
					//The ith particle moves a distance 'a' to the right and through boundary.
					x[i] += a;
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
		//fprintf(outdists, "%d ",t);
		for(i = 0; i < xL; i++){fprintf(outdists, "%d ", rho[i]);}
		fprintf(outdists,"\n");
		
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