/*
Non-particle-based diffusion simulation on 1D, homogenous lattice.
*/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define N 1.0
#define xC 15
#define xE 15
#define nU 11

int main(){
	//Index variables; for-loops and time-step limit.
	int i, j;
	int t, tmax = 20000;

	//Unit cell dimensions.
	int xU = xC + xE;

	//Total lattice dimensions.
	int xL = xU*nU;

	//Start lattice site (x,y); currently a cellular region.
	int xSP = (int)(xL/2.0) - 1;

	//Creating density distribution arrays.
	double rho_c[xL];
	double rho_n[xL];

	//Initializing density distribution arrays.
	for(i = 0; i < xL; i++){
		rho_c[i] = 0;
		rho_n[i] = 0;
	}
	//Setting density at start position of particles.
	rho_n[xSP] = N;

	//Stepping probabilities (arb. chosen) for intracellular.
	//Physical model: intracellular regions less diffusive (more viscous).
	double pnxi = 0.1, ppxi = 0.1, psxi = 1.0 - pnxi - ppxi;

	//For analytics.
	double sum_x, sum_x2, avg_x, avg_x2;
	int count = 0;

	//File naming and preparation. Why do we get a segmentation fault below?
	//char *path = "/home/paul/Documents/thesis/particle-diffusion/data/";
	//char *f1 = strcat(path,"TEST1.txt");
	//char *f2 = strcat(path,"TEST1-stats.txt");
	char *f1 = "/home/paul/Documents/thesis/particle-diffusion/data/H000.txt";
	char *f2 = "/home/paul/Documents/thesis/particle-diffusion/data/H000_stats.txt";
	FILE *outdists, *outstats;
	outdists = fopen(f1, "w");
	outstats = fopen(f2, "w");

	for(t = 0; t < tmax; t++){
		//Update the current time step density distribution (dd).

		//Resetting variables for analysis.
		sum_x = 0;
		sum_x2 = 0;
		avg_x = 0;
		avg_x2 = 0;
		count = 0; //to count number of lattice sites.

		for (i = 0; i < xL; i++){
			//Copy previous time step array into current array. Clear previous array (n).
			rho_c[i] = rho_n[i];
			rho_n[i] = 0;
		}
		
		//Now the current dd will be used to fill the next time-step array.

		for(i = 0; i < xL; i++){
			//Loop over all lattice sites and recalculate dds.
			
			// if(rho_c[i] == 0){
			// 	//No particles at current site, skip the remaining statements.
			// 	continue;
			// }

			if( (i != 0) & (i != xL - 1) ){
				//Not at an absolute boundary.

				rho_n[i] = psxi*rho_c[i] + 
				ppxi*rho_c[i-1] +
				pnxi*rho_c[i+1];
			}
			else{
				//At an absolute boundary.
				if(i == 0){
					//Zero probability to move -x (outside lattice).
					//If particle doesn't move +x, then in stays in place.
					rho_n[i] = (1.0-pnxi)*rho_c[i] +
					pnxi*rho_c[i+1];
				}
				else{
					//Zero probability to move +x (outside lattice).
					//If particle doesn't move -x, then in stays in place.
					rho_n[i] = (1.0-ppxi)*rho_c[i] +
					ppxi*rho_c[i-1];				
				}
			}
			//Note the difference from the MC computation for MSD.
			//Looping over lattice site here, not every particle.
			sum_x +=  (double)i*rho_n[i];
			sum_x2 += (double)i*(double)i*rho_n[i];
			count += 1;
		}

		for (i = 0; i < xL; i++){
			//Writing DD data to file.
			//Could maybe place this in the loop above.
			fprintf(outdists, "%f ", rho_n[i]);
			// //Calculating terms for MSD.
			// //This is done for every time step.
			// sum_x += i*rho_n[i];
			// sum_x2 = i*i*rho_n[i];
		}
		fprintf(outdists, "\n");

		//Writing mean-square-displacement data to file.
		avg_x = sum_x;
		avg_x2 = sum_x2;
		fprintf(outstats, "%f %f %f\n", avg_x, avg_x2, avg_x2 - avg_x*avg_x);

		//Print out a progress update.
		if(t%100 == 0){
		printf("%i\n",t);
		}
	}
}