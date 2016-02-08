/*
Diffusion simulation on 1D, inhomogenous (lattice)? Not particle based.

*/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define N 100	//can be set to whatever.
#define xC 15
#define xE 15
#define nU 3

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
	int ucolx; //unit cell of particle.
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
	double pie = 0.05;
	double pei = (pnxi/pnxe)*pie;

	//File naming and preparation. Why do we get a segmentation fault below?
	//char *path = "/home/paul/Documents/thesis/particle-diffusion/data/";
	//char *f1 = strcat(path,"TEST1.txt");
	//char *f2 = strcat(path,"TEST1-stats.txt");
	char *f1 = "/home/paul/Documents/thesis/particle-diffusion/data/t-3k_xU-3_pi-0.1_pe-0.2-pie-0.05.txt";
	char *f2 = "/home/paul/Documents/thesis/particle-diffusion/data/t-3k_xU-3_pi-0.1_pe-0.2-pie-0.05_stats.txt";
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

			//Some clarifying statements...
			ucolx = (int)(i/xU);
			modx = (int)(i - ucolx*xU);

			if( (i != 0) & (i != (xL - 1)) ){
				//Not at absolute boundary.
				
				if( modx < xC ){
					//Lattice site is in cellular region...
					if( (modx != 0) & (modx != (xC-1)) ){
						//Lattice site is not at cellular boundaries.
						rho_n[i - 1] += pnxi*rho_c[i];
						rho_n[i] += psxi*rho_c[i];
						rho_n[i + 1] += ppxi*rho_c[i];
					}
					else if(modx == 0){
						//Lattice site is at -x C-EC boundary.
						rho_n[i - 1] += pnxi*pie*rho_c[i]; //jump across boundary
						rho_n[i] += (1.0 - pnxi*pie - ppxi)*rho_c[i]; //stay in place
						rho_n[i + 1] += ppxi*rho_c[i]; //move right					
					}
					else{
						//Lattice site is at +x C-EC boundary.
						rho_n[i - 1] += pnxi*rho_c[i]; //move left
						rho_n[i] += (1.0 - pnxi - ppxi*pie)*rho_c[i]; //stay in place
						rho_n[i + 1] += ppxi*pie*rho_c[i]; //jump across boundary
					}
				}
				else{
					//Lattice site is in extracellular region...
					if( (modx != xC) & (modx != (xU-1)) ){
						//Lattice site is not at extracellular boundaries.
						rho_n[i - 1] += pnxe*rho_c[i];
						rho_n[i] += psxe*rho_c[i];
						rho_n[i + 1] += ppxe*rho_c[i];
					}
					else if(modx == xC){//maybe an i
						//Lattice site is at -x EC-C boundary.
						rho_n[i - 1] += pnxe*pei*rho_c[i]; //jump across boundary
						rho_n[i] += (1.0 - pnxe*pei - ppxe)*rho_c[i]; //stay in place
						rho_n[i + 1] += ppxe*rho_c[i];					
					}
					else{
						//Lattice site is at +x EC-C boundary.
						rho_n[i - 1] += pnxe*rho_c[i];
						rho_n[i] += (1.0 - pnxe - ppxe*pei)*rho_c[i]; //stay in place
						rho_n[i + 1] += ppxe*pei*rho_c[i]; //jump across boundary
					}
				}				
			}
			else{
				//At an absolute boundary.
				if(i == 0){
					//Zero probability to move -x (outside lattice).
					//If particle doesn't move +x, then in stays in place.
					rho_n[i] += (1.0 - ppxi - 0)*rho_c[i];
					rho_n[i + 1] += ppxi*rho_c[i]; 
				}
				else{
					//Zero probability to move +x (outside lattice).
					//If particle doesn't move -x, then in stays in place.
					rho_n[i - 1] += pnxi*rho_c[i];
					rho_n[i] += (1.0 - pnxi - 0)*rho_c[i];					
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