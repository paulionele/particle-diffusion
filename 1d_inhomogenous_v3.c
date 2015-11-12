/* 
Diffusion simulation on 1D, inhomogenous lattice. 

The code [should] work.
Boundary conditions, more specifically the probabilities of crossing the bounds 
are non-physical. The probablities of crossing from higher to lower diffusion
volumes cannot be the same. They are related by a formula, but not equal.

pei = something
pie = (De/Di)*pei = (pe/pi)*pei

The unit cell is built so that the overall cell lattice is non-symmetric. That's
something I need to fix.. I like symmetry.
*/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"

// CPP directives.
#define N 150000  //number of particles
#define xC 20    //x length of cellular space
#define xE 5  //x length of extracellular space
#define nU 3  //x number of unit cells

int main(){
	//Index variables for the for-loops and time-step limit.
	int i;
	int t, tmax = 5000;

	//x unit cell length (xU), total array length (xL) and centered start (SP).
	int xU = xC + xE;
	int xL = xU*nU;
	int SP = (int)(xL/2.0);

	///NEW VARIABLES DEFINED....
	int modx; //Modified x-pos. of ith particle. Returns: [0,xU-1]
	float testx = 0; //Possible future position of particle. Returns: [0,xL-1]
	float modtestx = 0;
	int state = 0;
	int teststate = 0;

	//Stepsize and probabilities for intracellular (0) and extracellular diffusion (1).
	//Probablities crossing from intra to extra (pie) or extra to intra (pei). 
	float a = 1.0;
	float pli = 0.3, pri = 0.3, psi = 1.0-pli-pri;
	float ple = 0.2, pre = 0.2, pse = 1.0-ple-pre;

	float pei = 0.05;
	float pie = (0.2/0.3)*pei; 

	//Random variable and variables for MSD calculations.
	float rnd;
	float sum_x, sum_x2, avg_x, avg_x2;

	//Shifted index and an array to hold particle positions.
	int sxi; //is this still sxi? Does that name make sense?
	int ucop; //unit cell of particle.
	int ucot; //unit cell of test.
	float x[N] = {[0 ... (N-1)] = SP}; //designated initializer to set start positions
	//for(i = 0; i < N; i++){x[i] = SP;} //initialize with a rolling loop instead.


	//File naming and preparation. Why do we get a segmentation fault below?
	//char *path = "/home/paul/Documents/thesis/particle-diffusion/data/";
	//char *f1 = strcat(path,"TEST1.txt");
	//char *f2 = strcat(path,"TEST1-stats.txt");
	char *f1 = "/home/paul/Documents/thesis/particle-diffusion/data/TEST1.txt";
	char *f2 = "/home/paul/Documents/thesis/particle-diffusion/data/TEST1-stats.txt";
	FILE *outdists, *outstats;
	outdists = fopen(f1, "w");
	outstats = fopen(f2, "w");
	
	//Need to have rho pre-defined for use in the main loop.
	int *rho;
	rho = malloc(xL*sizeof(int));
	
	//No need to touch: create density array and print out first distribution.
	for(i = 0; i < N; i++){sxi = (int)x[i]; rho[sxi]++;}
	for(i = 0; i < xL; i++){fprintf(outdists, "%d ", rho[i]);}
	fprintf(outdists, "\n");

	//--------------------------------------------------------------------------
	//Below is where the main work is done.
	for (t = 0; t < tmax; t++){
		sum_x = 0;
		sum_x2 = 0;

		//Set density array elements to zero.
		for (i = 0; i < xL; i++){
			rho[i] = 0;
		}

		for(i = 0; i < N; i++){
			//Loop over all particles.
			//x[i]: [0, n*xU-1 = xL -1]
			//ucop: [0,n-1]
			//modx: [0,xU-1]
			ucop = (int)(x[i]/xU);
			modx = x[i] - ucop*xU;
			if(modx < xC){
				//Particle is currently in intracellular region.
				state = 1;
			}
			else{
				//Particle is currently in extracellular region.
				state = 0;
			}

			if(state == 1){
				//Apply intracellular conditions to test particle.
				rnd = (float)rand()/(float)RAND_MAX;
				if (rnd < pli){
					//Generate position of particle IF to move left.
					//Off-lattice if testx < 0, x[i] = 0 for that to occur.
					testx = x[i] - a;
				}
				else if(rnd < (pli + pri)){
					//Generate position of partcle IF to move right.
					testx = x[i] + a;
				}
				else{
					//Generate position of particle if to stay.
					testx = x[i];
				}
			}
			else{
				//Apply extracellular conditions to test particle.
				rnd = (float)rand()/(float)RAND_MAX;
				if (rnd < ple){
					//Generate position of particle if to move left.
					testx = x[i] - a;
				}
				else if(rnd < (ple + pre)){
					//Generate position of particle if to move right.
					//Off-lattice if testx > (xL - 1)
					testx = x[i] + a;
				}
				else{
					testx = x[i];
				}
			}

			if(testx >= 0 && testx < xL && testx != x[i]){
				//Det. motion of the particle depending on test position.
				//If at absolute lattice limits, do not do anything.
				//ucot: [0,n-1]
				//modtestx: [0,xU-1]
				ucot = (int)(testx/xU);
				modtestx = testx - ucot*xU;
				if(modtestx < xC){
					teststate = 1;
				} 
				else{
					teststate = 0;
				}
				
				//Has the particle moved to a different cell region?
				if(state == teststate){
					//Particle has not moved to different cell region.
					x[i] = testx;
				}
				else if(state == 1 && teststate == 0){
					//Particle is in intracellular region and at boundary.
					//Draw rnd to determine movement.
					rnd = (float)rand()/(float)RAND_MAX;
					if(rnd < pie)
						//Particle will cross boundary into extracellular region.
						x[i] = testx;
				}
				else{
					//Particle is in extracellular region and at boundary.
					//Draw rnd to determine movement.
					rnd = (float)rand()/(float)RAND_MAX;
					if(rnd < pei)
						//Particle will cross boundary in intracellular region.
						x[i] = testx;					
				}
			}

			//No need to touch anything below this line.
			sum_x += x[i];
			sum_x2 += x[i]*x[i];

			/*Shift the particle position and increment particle count at that
			position rho[sxi] by 1 unit.*/
			sxi = (int)(x[i]);
			rho[sxi]++; 
		}
		
		//Writing cell density data to file. First column is time.
		//fprintf(outdists, "%d ",t);
		for(i = 0; i < xL; i++){fprintf(outdists, "%d ", rho[i]);}
		fprintf(outdists,"\n");
		
		//Writing mean-square-displacement data to file.
		avg_x = sum_x/(float)N;
		avg_x2 = sum_x2/(float)N;
		fprintf(outstats, "%f %f %f\n", avg_x, avg_x2, avg_x2-avg_x*avg_x);
	}
	//return 0;
}

// get_argv(){
// 	//Get command-line arguements.
// 	//
// }