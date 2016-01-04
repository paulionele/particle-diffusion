/* 
Diffusion simulation on 1D, inhomogenous lattice. 

The code [should] work.
Boundary conditions, more specifically the probabilities of crossing the bounds 
are non-physical. The probablities of crossing from higher to lower diffusion
volumes cannot be the same. They are related by a formula, but not equal.

Smaller probability of transition from more diffusive (e) to less diffusive (i).
pei = (Di/De)*pie

The unit cell is built so that the overall cell lattice is non-symmetric.
*/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define N 500000 //number of particles
#define xC 11    //x length of cellular space
#define xE 11  //x length of extracellular space
#define nU 13  //x number of unit cells

int main(){
	//Index variables for the for-loops and time-step limit.
	int i;
	int t, tmax = 10000;

	//x unit cell length (xU), total array length (xL) and centered start (SP).
	int xU = xC + xE;
	int xL = xU*nU;
	int SP = (int)(xL/2.0) - (int)(xC/2.0) - 1;

	int modx; //Modified x-pos. of ith particle. Returns: [0,xU-1]
	double testx = 0; //Possible future state of particle. Returns: [0,xL-1]
	double modtestx = 0;
	int state = 0;
	int teststate = 0;

	int sxi; //used in converting ith particle position to integer.
	int ucop; //unit cell of particle.
	int ucot; //unit cell of test.
	double x[N] = {[0 ... (N-1)] = SP}; //designated initializer to set start positions
	//for(i = 0; i < N; i++){x[i] = SP;} //initialize with a rolling loop instead.

	//Particle density distribution array.
	int rho[xL];

	//Stepsize and stepping probablities (arb. chosen); equal for unbiased motion.
	double a = 1.0;
	double pli = 0.2, pri = 0.2, psi = 1.0-pli-pri; //Intracellular less diffusive.
	double ple = 0.3, pre = 0.3, pse = 1.0-ple-pre; //Extracellular more diffusive.
	//Coupled probablities crossing from intra to extra (pie) or extra to intra (pei). 
	double pie = 0.05; //Arbitrarily chosen.
	double pei = (pli/ple)*pie; 

	//Random variable and variables for MSD calculations.
	double rnd;
	double sum_x, sum_x2, avg_x, avg_x2;

	//File naming and preparation. Why do we get a segmentation fault below?
	//char *path = "/home/paul/Documents/thesis/particle-diffusion/data/";
	//char *f1 = strcat(path,"TEST1.txt");
	//char *f2 = strcat(path,"TEST1-stats.txt");
	char *f1 = "/home/paul/Documents/thesis/particle-diffusion/data/t000.txt";
	char *f2 = "/home/paul/Documents/thesis/particle-diffusion/data/t000_stats.txt";
	FILE *outdists, *outstats;
	outdists = fopen(f1, "w");
	outstats = fopen(f2, "w");
	
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
				rnd = (double)rand()/(double)RAND_MAX;
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
				rnd = (double)rand()/(double)RAND_MAX;
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
					rnd = (double)rand()/(double)RAND_MAX;
					if(rnd < pie)
						//Particle will cross boundary into extracellular region.
						x[i] = testx;
				}
				else{
					//Particle is in extracellular region and at boundary.
					//Draw rnd to determine movement.
					rnd = (double)rand()/(double)RAND_MAX;
					if(rnd < pei)
						//Particle will cross boundary in intracellular region.
						x[i] = testx;					
				}
			}

			//No need to touch anything below this line.
			sum_x += x[i];
			sum_x2 += x[i]*x[i];

			/*Convert particle position to integer and increment particle count 
			at that	position rho[sxi] by 1 unit.*/
			sxi = (int)(x[i]);
			rho[sxi]++; 
		}
		
		//Writing density distribution data to file.
		for(i = 0; i < xL; i++){fprintf(outdists, "%d ", rho[i]);}
		fprintf(outdists,"\n");
		
		//Writing mean-square-displacement data to file.
		avg_x = sum_x/(double)N;
		avg_x2 = sum_x2/(double)N;
		fprintf(outstats, "%f %f %f\n", avg_x, avg_x2, avg_x2-avg_x*avg_x);
	}
	//return 0;
}

// get_argv(){
// 	//Get command-line arguements.
// 	//
// }