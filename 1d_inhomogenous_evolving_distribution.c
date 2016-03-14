/*
Non-particle-based diffusion simulation on 1D, inhomogenous lattice.
Probabilistic method of advancing the distribution not the same as
the 1D homogenous program or 2D heterogenous.
*/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define N 1.0	//can be set to whatever.
#define xC 15
#define xE 15
#define nU 11

int main(){
	//Index variables; for-loops and time-step limit.
	int i, j;
	int t, tmax = 20000;

	//Unit cell and absolute dimensions.
	int xU = xC + xE;
	int xL = xU*nU;

	//Start lattice site (x); currently a cellular region.
	int xSP = (int)(xL/2.0) - (int)(xC/2.0) - 1;
	//int xSP = 1;

	//Variables for determining movement behavior.
	int ucolx; //unit cell of lattice.
	int modx; //modified x-position.
	//Used in converting position to integer.
	int sxi;

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
	double pnxe = 0.2, ppxe = 0.2, psxe = 1.0 - pnxe - ppxe;
	//Transition probabilities between cell types.
	double pie = 0.025;
	double pei = (pnxi/pnxe)*pie;

	//For analytics.
	double sum_x, sum_x2, avg_x, avg_x2;
	int count = 0;

	//File naming and preparation. Why do we get a segmentation fault below?
	//char *path = "/home/paul/Documents/thesis/particle-diffusion/data/";
	//char *f1 = strcat(path,"TEST1.txt");
	//char *f2 = strcat(path,"TEST1-stats.txt");
	char *f1 = "/home/paul/Documents/thesis/particle-diffusion/data/IH000.txt";
	char *f2 = "/home/paul/Documents/thesis/particle-diffusion/data/IH000_stats.txt";
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

			//Determining unit cell of current lattice site and position within.
			ucolx = (int)(i/xU);
			modx = (int)(i - ucolx*xU);

			if( (i != 0) & (i != (xL - 1)) ){
				//Not at absolute boundary.
				
				if( modx < xC ){
					//Lattice site is in cellular region...
					if( (modx != 0) & (modx != (xC-1)) ){
						//Lattice site is not at cellular boundaries.
						rho_n[i] = psxi*rho_c[i] +
						ppxi*rho_c[i-1] +
						pnxi*rho_c[i+1];
					}
					else if(modx == 0){
						//Lattice site is at -x C-EC boundary.
						rho_n[i] = (1.0-pie*pnxi-ppxi)*rho_c[i] +
						pei*ppxe*rho_c[i-1] +
						pnxi*rho_c[i+1];				
					}
					else{
						//Lattice site is at +x C-EC boundary.
						rho_n[i] = (1.0-pie*ppxi-pnxi)*rho_c[i] +
						ppxi*rho_c[i-1] +
						pei*pnxe*rho_c[i+1];
					}
				}
				else{
					//Lattice site is in extracellular region...
					if( (modx != xC) & (modx != (xU-1)) ){
						//Lattice site is not at extracellular boundaries.
						rho_n[i] = psxe*rho_c[i] +
						ppxe*rho_c[i-1] +
						pnxe*rho_c[i+1];
					}
					else if(modx == xC){
						//Lattice site is at -x EC-C boundary.
						rho_n[i] = (1.0-pei*pnxe-ppxe)*rho_c[i] +
						pie*ppxi*rho_c[i-1] +
						pnxe*rho_c[i+1];					
					}
					else{
						//Lattice site is at +x EC-C boundary.
						rho_n[i] = (1.0-pei*ppxe-pnxe)*rho_c[i] +
						ppxe*rho_c[i-1] +
						pie*pnxi*rho_c[i+1];
					}
				}				
			}
			else{
				//At an absolute boundary.
				if(i == 0){
					//Zero probability to move -x (outside lattice).
					//If particle doesn't move +x, then in stays in place.
					rho_n[i] = (1.0 - ppxi - 0)*rho_c[i] +
					pnxi*rho_c[i+1];
				}
				else{
					//Zero probability to move +x (outside lattice).
					//If particle doesn't move -x, then in stays in place.
					rho_n[i] = (1.0 - pnxe - 0)*rho_c[i] +
					ppxe*rho_c[i-1];				
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