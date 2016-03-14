/*
Numerically solving the diffusion equation, on a lattice,
boundary conditions applied.
*/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define N 1.0	//can be set to whatever.
#define xC 15
#define xE 15
#define yC 15
#define yE 15
#define nU 11
#define mU 1

int main(){
	//Index variables; for-loops and time-step limit.
	int i, j;
	int mi, mj;
	int t, tmax = 20000;

	//Unit cell dimensions.
	int xU = xC + xE;
	int yU = yC + yE;
	//Total lattice/array dimensions.
	int xL = xU*nU;
	int yL = yU*mU;

	//Start lattice site (x,y); currently a cellular region.
	int xSP = (int)(xL/2.0) - (int)(xC/2.0) - 1;
	int ySP = (int)(yL/2.0) - (int)(yC/2.0) - 1;
	
	//Unit cell of lattice and modified position of lattice.
	int ucolx;
	int ucoly;
	int modx;
	int mody;

	//Used in converting position to integer.
	int sxi;
	int syi;

	//Current time step and next time step density distribution arrays.
	double rho_c[yL+1][xL+1];
	double rho_n[yL+1][xL+1];
	double mask[yL+1][xL+1];
	
	//Probability distribution arrays.
	double ppx[yL+1][xL+1];
	double pnx[yL+1][xL+1];
	double ppy[yL+1][xL+1];
	double pny[yL+1][xL+1];
	double ps[yL+1][xL+1];

	//Initializing density distribution arrays.
	for (j = 0; j <= yL; j++){
		for(i = 0; i <= xL; i++){
			rho_c[j][i] = 0;
			rho_n[j][i] = 0;
			ppx[j][i]=0;
			pnx[j][i]=0;
			ppy[j][i]=0;
			pny[j][i]=0;
			ps[j][i]=0;
		}
	}

	//Setting density at start position of particles.
	rho_n[ySP][xSP] = N;

	//Stepping probabilities (arb. chosen) for intracellular.
	//Physical model: intracellular regions less diffusive (more viscous).
	double pnxi = 0.05, ppxi = 0.05, psxi = 1.0 - pnxi - ppxi;
	double pnyi = 0.05, ppyi = 0.05, psyi = 1.0 - pnyi - ppyi;
	//Stepping probabilities (arb. chosen) for extracellular.
	//Physical model: extracellular regions more diffusive (less viscous).
	double pnxe = 0.2, ppxe = 0.2, psxe = 1.0 - pnxe - ppxe;
	double pnye = 0.2, ppye = 0.2, psye = 1.0 - pnye - ppye;
	//Transition probabilities between cell types.
	double pie = 0.01;
	double pei = (pnxi/pnxe)*pie;

	//For analytics.
	double sum_x, sum_x2, avg_x, avg_x2;
	double sum_y, sum_y2, avg_y, avg_y2;
	int count = 0;

	//File naming and preparation. Why do we get a segmentation fault below?
	//char *path = "/home/paul/Documents/thesis/particle-diffusion/data/";
	//char *f1 = strcat(path,"TEST1.txt");
	//char *f2 = strcat(path,"TEST1-stats.txt");
	char *f1 = "/home/paul/Documents/thesis/particle-diffusion/2D/2D-data/0FD_t-30k_N-1_nU-21_pi-0.05_pe-0.2_pie-0.01.txt";
	char *f2 = "/home/paul/Documents/thesis/particle-diffusion/2D/2D-data/0FD_t-600k_N-1_nU-41_pi-0.05_pe-0.2_pie-0.01_stats.txt";
	FILE *outdists, *outstats;
	outdists = fopen(f1, "w+");
	outstats = fopen(f2, "w+");

	//Create mask of lattice point types.
	for(j = 0; j <= yL; j++){
		for(i = 0; i <= xL; i++){

			//Set entire system to type-2.
			mask[j][i] = 2;
			
			ucolx = (int)(i/xU);
			ucoly = (int)(j/yU);
			modx = (int)(i - ucolx*xU);
			mody = (int)(j - ucoly*yU);
			
			//Find and set those intracellular sites specifically to type-1.
			if( (modx < xC) && (mody < yC) ){
				mask[j][i] = 1;	
			}
		}
	}
	
	//Set the border sites to type-0. No lattice points used here.
	//j = 0;
	for(i = 0; i <= xL; i++){
		mask[0][i] = 0;
	}
	//j=yL;
	for(i = 0; i <= xL; i++){
		mask[yL][i] = 0;
	}	
	//i=0;
	for(j = 0; j <= yL; j++){
		mask[j][0] = 0;
	}
	//i=xL;
	for(j = 0; j <= yL; j++){
		mask[j][xL] = 0;
	}

	//Creating the probability matrix.
	for(j = 0; j <= yL; j++){
		for(i = 0; i <= xL; i++){
			//Loop over all the lattice points.

			if(mask[j][i] != 0){
				//Lattice site is not a border lattice.

				if(mask[j][i] == 1){
					//Intracellular lattice point.
					//mi and mj are modified indicies.

					mi = i-1;
					mj = j;
					if(mask[mj][mi] == 0){
						//No particles move -x.
						pnx[j][i] = 0;
					}
					else if(mask[mj][mi] == 1){
						//If adjacent lattice of same type, probability to move left.
						pnx[j][i] = pnxi;
					}
					else if(mask[mj][mi] == 2){
						//If adjacent lattice of diff type, prod prob to move right.
						pnx[j][i] = pnxi*pie;
					}
					else{
						printf("Error Mask at i-1!\n");
					}
					
					mi = i+1;
					mj = j;
					if(mask[mj][mi] == 0){
						ppx[j][i] = 0;
					}
					else if(mask[mj][mi] == 1){
						ppx[j][i] = ppxi;
					}
					else if(mask[mj][mi] == 2){
						ppx[j][i] = ppxi*pie;
					}
					else{
						printf("Error Mask at i+1!\n");
					}

					mi = i;
					mj = j-1;
					if(mask[mj][mi] == 0){
						pny[j][i] = 0;
					}
					else if(mask[mj][mi] == 1){
						pny[j][i] = pnyi;
					}
					else if(mask[mj][mi] == 2){
						pny[j][i] = pnyi*pie;
					}
					else{
						printf("Error Mask at j-1!\n");
					}
					
					mi = i;
					mj = j+1;
					if(mask[mj][mi] == 0){
						ppy[j][i] = 0;
					}
					else if(mask[mj][mi] == 1){
						ppy[j][i] = ppyi;
					}
					else if(mask[mj][mi] == 2){
						ppy[j][i] = ppyi*pie;
					}
					else{
						printf("Error Mask at j+1!\n");
					}

					//The probability of staying at current lattice site.
					//Subtraction of all the probabilities of motion at current site.
					ps[j][i]=1.0-pnx[j][i]-ppx[j][i]-pny[j][i]-ppy[j][i];					

				}
				else if(mask[j][i] == 2){
					//Extracellular lattice points.
					//mi and mj are modified indicies.
					
					mi = i-1;
					mj = j;
					if(mask[mj][mi] == 0){
						pnx[j][i] = 0;
					}
					else if(mask[mj][mi] == 1){
						pnx[j][i] = pnxe*pei;
					}
					else if(mask[mj][mi] == 2){
						pnx[j][i] = pnxe;
					}
					else{
						printf("Error Mask at i-1!\n");
					}
					
					mi = i+1;
					mj = j;
					if(mask[mj][mi] == 0){
						ppx[j][i] = 0;
					}
					else if(mask[mj][mi] == 1){
						ppx[j][i] = ppxe*pei;
					}
					else if(mask[mj][mi] == 2){
						ppx[j][i] = ppxe;
					}
					else{
						printf("Error Mask at i+1!\n");
					}

					mi = i;
					mj = j-1;
					if(mask[mj][mi] == 0){
						pny[j][i] = 0;
					}
					else if(mask[mj][mi] == 1){
						pny[j][i] = pnye*pei;
					}
					else if(mask[mj][mi] == 2){
						pny[j][i] = pnye;
					}
					else{
						printf("Error Mask at j-1!\n");
					}
					
					mi = i;
					mj = j+1;
					if(mask[mj][mi] == 0){
						ppy[j][i] = 0;
					}
					else if(mask[mj][mi] == 1){
						ppy[j][i] = ppye*pei;
					}
					else if(mask[mj][mi] == 2){
						ppy[j][i] = ppye;
					}
					else{
						printf("Error Mask at j+1!\n");
					}
					
					//The probability of staying at the current lattice site.
					ps[j][i]=1.0-pnx[j][i]-ppx[j][i]-pny[j][i]-ppy[j][i];
				}
				else{
					printf("Error Mask!\n");
				}
			}
		}
	}
	
	//The above is only calculated once.
	for(t = 0; t < tmax; t++){
		//Update the current time step density distribution (dd).
		//Now the current dd will be used to fill the next time-step array.

		//Resetting variables for analysis.
		sum_x = 0;
		sum_x2 = 0;
		avg_x = 0;
		avg_x2 = 0;
		count = 0; //will equal the number of (real?) lattice site.

		//Updating density distribution.
		for (j = 0; j <= yL; j++){
			for(i = 0; i <= xL; i++){
				rho_c[j][i] = rho_n[j][i];
				rho_n[j][i] = 0;
			}
		}

		//Loop over all lattice sites and recalculate dds.
		float rhot = 0;
		for(j = 0; j <= yL; j++){
			for(i = 0; i <= xL; i++){
				if(mask[j][i] != 0){
					rho_n[j][i] = 
						ps[j][i]*rho_c[j][i] +
						ppx[j][i-1]*rho_c[j][i-1] +
						pnx[j][i+1]*rho_c[j][i+1] +
						ppy[j-1][i]*rho_c[j-1][i] +
						pny[j+1][i]*rho_c[j+1][i];

					//Note the difference from the MC computation for MSD.
					//Looping over lattice site here, not every particle.
					//Also, don't include the border lattice points.
					//pow((double)(i), 2)
					sum_x +=  (double)i*rho_n[j][i];
					sum_x2 += (double)i*(double)i*rho_n[j][i];
					count += 1;
				}

				//Total after sum should equal ~1.0 (if N = 1.0).
				rhot += rho_n[j][i];
				
				// if(rho_n[j][i]>1.0){
				// 	//printf("Arg!          %i %i %f\n",i,j,rho_n[j][i]);
				// 	//getchar();
				// }
			}
		}

		//Print out a progress update.
		if(t%100 == 0){
			printf("%i %f\n",t,rhot);
		}
		
		//Writing dd data to file.
		for (j = 0; j <= yL; j++){
			for(i = 0; i <= xL; i++){
				fprintf(outdists, "%f ", rho_n[j][i]);
				//fprintf(outdists, "%f ", mask[j][i]);
			}
			fprintf(outdists, "\n");
		}
		fprintf(outdists, "\n");

		//Writing MSD data to file, done for every time step.
		//Don't need to divide by count (the number of lattice sites)
		// unless N != 1.
		avg_x = sum_x;
		avg_x2 = sum_x2;
		fprintf(outstats, "%lf %lf %lf\n", avg_x, avg_x2, avg_x2-avg_x*avg_x);
	}
}