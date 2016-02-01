/*
Diffusion simulation on 1D, homogenous (lattice)? Not particle based.

*/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define N 100
#define xC 11
#define xE 11
#define nU 3

int main(){
	//Index variables; for-loops and time-step limit.
	int i, j;
	int t, tmax = 3000;

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

	//File naming and preparation. Why do we get a segmentation fault below?
	//char *path = "/home/paul/Documents/thesis/particle-diffusion/data/";
	//char *f1 = strcat(path,"TEST1.txt");
	//char *f2 = strcat(path,"TEST1-stats.txt");
	char *f1 = "/home/paul/Documents/thesis/particle-diffusion/data/t-3k_xU-3_p-0.1.txt";
	char *f2 = "/home/paul/Documents/thesis/particle-diffusion/data/t-3k_xU-3_p-0.1_stats.txt";
	FILE *outdists, *outstats;
	outdists = fopen(f1, "w");
	outstats = fopen(f2, "w");

	for(t = 0; t < tmax; t++){
		//Update the current time step density distribution (dd).

		for (i = 0; i < xL; i++){
			//Copy previous time step array into current array. Clear previous array (n).
			rho_c[i] = rho_n[i];
			rho_n[i] = 0;
		}
		
		//Now the current dd will be used to fill the next time-step array.

		for(i = 0; i < xL; i++){
			//Loop over all lattice sites and recalculate dds.
			
			if(rho_c[i] == 0){
				//No particles at current site, skip the remaining statements.
				continue;
			}

			if( (i != 0) & (i != xL - 1) ){
				//Not at boundary. Advance the density distribution.

				rho_n[i - 1] += pnxi*rho_c[i];
				rho_n[i] += psxi*rho_c[i];
				rho_n[i + 1] += ppxi*rho_c[i]; 
			}
			else{
				//At a boundary. What behaviour should we have?
				if(i == 0){
					//Half stay in place, half move away???
					rho_n[i] += 0.5*rho_c[i];
					rho_n[i + 1] += 0.5*rho_c[i]; 
				}
				else{
					//Half stay in place, half move away???
					rho_n[i - 1] += 0.5*rho_c[i];
					rho_n[i] += 0.5*rho_c[i];					
				}
			}
		}
		//Writing dd data to file.
		for (i = 0; i < xL; i++){
			fprintf(outdists, "%f ", rho_n[i]);
		}
		fprintf(outdists, "\n");
	}
}