/*
Numerically solving the diffusion equation, on a lattice,
boundary conditions applied.

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
	double rho_c[yL][xL];
	double rho_n[yL][xL];

	//Initializing density distribution arrays.
	for (j = 0; j < yL; j++){
			for(i = 0; i < xL; i++){
				rho_c[j][i] = 0;
				rho_n[j][i] = 0;
			}
		}

	//Setting density at start position of particles.
	rho_n[ySP][xSP] = N;

	//Stepping probabilities (arb. chosen) for intracellular.
	//Physical model: intracellular regions less diffusive (more viscous).
	double pnxi = 0.2, ppxi = 0.2, psxi = 1.0 - pnxi - ppxi;
	double pnyi = 0.2, ppyi = 0.2; psyi = 1.0 - pnyi - ppyi;
	//Stepping probabilities (arb. chosen) for extracellular.
	//Physical model: extracellular regions more diffusive (less viscous).
	double pnxe = 0.4, ppxe = 0.4; psxe = 1.0 - pnxe - ppxe;
	double pnye = 0.4, ppye = 0.4; psye = 1.0 - pnye - ppye;
	//Transition probabilities between cell types.
	double pie = 0.025;
	double pei = (pnxi/pnxe)*pie;

	//File naming and preparation. Why do we get a segmentation fault below?
	//char *path = "/home/paul/Documents/thesis/particle-diffusion/data/";
	//char *f1 = strcat(path,"TEST1.txt");
	//char *f2 = strcat(path,"TEST1-stats.txt");
	char *f1 = "/home/paul/Documents/thesis/particle-diffusion/2D/2D-data/t-3k_xU-3_pi-0.1_pe-0.2-pie-0.05.txt";
	char *f2 = "/home/paul/Documents/thesis/particle-diffusion/2D/2D-data/t-3k_xU-3_pi-0.1_pe-0.2-pie-0.05_stats.txt";
	FILE *outdists, *outstats;
	outdists = fopen(f1, "w");
	outstats = fopen(f2, "w");

	for(t = 0; t < tmax; t++){
		
		//Update the current time step density distribution (dd).
		//Now the current dd will be used to fill the next time-step array.
		for (j = 0; j < yL; j++){
			for(i = 0; i < xL; i++){
				rho_c[j][i] = rho_n[j][i];
				rho_n[j][i] = 0;
			}
		}

		//Loop over all lattice sites and recalculate dds.
		for(j = 0; j < yL; j++){
			
			for(i = 0; i < xL; i++){
				
				if(rho_c[j][i] == 0){
					//No particles at current site, skip the remaining statements.
					//Probably doesn't provide much speed improvement in long run?!
					continue;
				}

				//Some clarifying statements to be added...
				ucolx = (int)(i/xU);
				ucoly = (int)(j/yU)
				modx = (int)(i - ucolx*xU);
				mody = (int)(j - ucoly*yU);

				if( ((i != 0) & (i != (xL - 1))) & ((j != 0) & (j != (yL - 1)))){
					//Not at absolute boundary.
					
					if(mody < yC){
						/*
						Lattice site is in top half of 2D array and may be in either an
						intracelllar or extracellular region.
						*/
//-------------
						if( modx < xC ){
							//Lattice site is in cellular region.
							//Testing for 3 different boundaries and 2 corner sites.
							if( ((modx != 0) & (modx != (xC-1))) & ((mody != 0) & (mody != (yC-1))) ){
								//Lattice site is not at cellular boundaries.
								rho_n[j][i] = 
								(1.0-ppxi*psyi-pnxi*psyi-psxi*ppyi-psxi*pnyi-ppxi*ppyi-pnxi*pnyi-ppxi*pnyi-pnxi*ppyi)*rho_c[j][i] + 
								ppxi*psyi*rho_c[j][i-1] +
								pnxi*psyi*rho_c[j][i+1] +
								psxi*ppyi*rho_c[j-1][i] +
								psxi*pnyi*rho_c[j+1][i] +
								ppxi*ppyi*rho_c[j-1][i-1] +
								pnxi*pnyi*rho_c[j+1][i+1] +
								ppxi*pnyi*rho_c[j+1][i-1] +
								pnxi*ppyi*rho_c[j-1][i+1]
							}
							else if( (mody == 0) & (modx != 0) & (modx != (xC-1)) ){
								//At the -y boundary, not at corner.
								//This case should not execute here, since only 1 unit cell in y.
								printf("Error, -y boundary.\n");
							}
							else if( (mody == (yC-1) & (modx != 0) & (modx != (xC-1)) ){
								//At the +y boundary, not at corner.
								rho_n[j][i] = 
								(1.0-ppxi*psyi-pnxi*psyi-psxi*ppyi-pei*psxe*pnye-ppxi*ppyi-pei*pnxe*pnye-pei*ppxe*pnye-pnxi*ppyi)*rho_c[j][i] + 
								ppxi*psyi*rho_c[j][i-1] +
								pnxi*psyi*rho_c[j][i+1] +
								psxi*ppyi*rho_c[j-1][i] +
								pei*psxe*pnye*rho_c[j+1][i] +
								ppxi*ppyi*rho_c[j-1][i-1] +
								pei*pnxe*pnye*rho_c[j+1][i+1] +
								pei*ppxe*pnye*rho_c[j+1][i-1] +
								pnxi*ppyi*rho_c[j-1][i+1]
							}
							else if( (modx == 0) & (mody != 0) & (mody != (yC-1)) ){
								//At the -x boundary, not at corner.
								rho_n[j][i] = 
								(1.0-pei*ppxe*psye-pnxi*psyi-psxi*ppyi-psxi*pnyi-pei*ppxe*ppye-pnxi*pnyi-pei*ppxe*pnye-pnxi*ppyi)*rho_c[j][i] + 
								pei*ppxe*psye*rho_c[j][i-1] +
								pnxi*psyi*rho_c[j][i+1] +
								psxi*ppyi*rho_c[j-1][i] +
								psxi*pnyi*rho_c[j+1][i] +
								pei*ppxe*ppye*rho_c[j-1][i-1] +
								pnxi*pnyi*rho_c[j+1][i+1] +
								pei*ppxe*pnye*rho_c[j+1][i-1] +
								pnxi*ppyi*rho_c[j-1][i+1]
							}
							else if( (modx == (xC-1)) & (mody != 0) & (mody != (yC-1)) ){
								//At the +x boundary, not at a corner.
								rho_n[j][i] = 
								(1.0-ppxi*psyi-pei*pnxe*psye-psxi*ppyi-psxi*pnyi-ppxi*ppyi-pei*pnxe*pnye-ppxi*pnyi-pei*pnxe*ppye)*rho_c[j][i] + 
								ppxi*psyi*rho_c[j][i-1] +
								pei*pnxe*psye*rho_c[j][i+1] +
								psxi*ppyi*rho_c[j-1][i] +
								psxi*pnyi*rho_c[j+1][i] +
								ppxi*ppyi*rho_c[j-1][i-1] +
								pei*pnxe*pnye*rho_c[j+1][i+1] +
								ppxi*pnyi*rho_c[j+1][i-1] +
								pei*pnxe*ppye*rho_c[j-1][i+1]
							}
							else if( modx == 0 & mody == 0){
								//At the -x,-y corner.
								printf("Error, -x,-y corner.\n");
							}
							else if( modx == (xC-1) & mody == 0){
								//At the +x,-y corner.
								printf("Error, +x,-y corner.\n");
							}
							else if( modx == (xC-1) & mody == (yC-1) ){
								//At the +x,+y corner.
								rho_n[j][i] = 
								(1.0-ppxi*psyi-pei*pnxe*psye-psxi*ppyi-pei*psxe*pnye-ppxi*ppyi-pei*pnxe*pnye-pei*ppxe*pnye-pei*pnxe*ppye)*rho_c[j][i] + 
								ppxi*psyi*rho_c[j][i-1] +
								pei*pnxe*psye*rho_c[j][i+1] +
								psxi*ppyi*rho_c[j-1][i] +
								pei*psxe*pnye*rho_c[j+1][i] +
								ppxi*ppyi*rho_c[j-1][i-1] +
								pei*pnxe*pnye*rho_c[j+1][i+1] +
								pei*ppxe*pnye*rho_c[j+1][i-1] +
								pei*pnxe*ppye*rho_c[j-1][i+1]
							}
							else if( modx == 0 & mody == (yC-1) ){
								//At the -x,+y corner.
								rho_n[j][i] = 
								(1.0-pei*ppxe*psye-pnxi*psyi-psxi*ppyi-pei*psxe*pnye-pei*ppxe*ppye-pei*pnxe*pnye-pei*ppxe*pnye-pnxi*ppyi)*rho_c[j][i] + 
								pei*ppxe*psye*rho_c[j][i-1] +
								pnxi*psyi*rho_c[j][i+1] +
								psxi*ppyi*rho_c[j-1][i] +
								pei*psxe*pnye*rho_c[j+1][i] +
								pei*ppxe*ppye*rho_c[j-1][i-1] +
								pei*pnxe*pnye*rho_c[j+1][i+1] +
								pei*ppxe*pnye*rho_c[j+1][i-1] +
								pnxi*ppyi*rho_c[j-1][i+1]
							}
							else{
								printf("Lattice site location error in cellular region.\n")
							}
						}
//---------------
						else{
							//Lattice site is in upper half extracellular region.
							//Testing for 2 different boundaries and 2 corner sites.

						}
					}
					else{
						/*
						The lattice site is in the lower half of the 2D array and is only
						in extracellular region. Refer to cell model.
						*/
					}
				}
				else{
					//At an absolute boundary.					
				}
			} //no touch
		} //no touch


		//Writing dd data to file.
		for (i = 0; i < xL; i++){
			fprintf(outdists, "%f ", rho_n[i]);
		}
		fprintf(outdists, "\n");
	}
}