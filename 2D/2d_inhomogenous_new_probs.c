/* 
DO NOT USE: This program contains the old particle movement algorithm. Future
modifications should be done on 2d_inhmogenous.c

Diffusion simulation on 2D, inhomogeneous lattice. 

Representation of lattice inhomogeneity.
101010
000000

NOTE: Due to the way determination of particle states is handled, usage of
more than one unit cell in column (mU > 1) may produce unexpected results. This
issue is derived from the fact that in the y-direction, there are two kinds of
unit cells. The current handling of the problem was intended only as a quick
solution to get the program running producing data for modeling.
*/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

// Parameters for characterization of entire 2D inhomogeneous lattice.
#define N 500000	//number of particles
#define xC 11		//x length of cellular space
#define yC 11		//y length of cellular space
#define xE 11		//x length of extracellular space
#define yE 11		//y length of extracellular space
#define nU 3		//number of unit cells in row
#define mU 1        //number of unit cells in column
#define mR 2		//number of rows (NOT USED)
#define nR 0        //number of columns (NOT USED)

int main(){
	//Index variables; for-loops and time-step limit.
	int i, j, n;
	int t, tmax = 10000;

	//Unit cell dimensions.
	int xU = xC + xE;
	int yU = yC + yE;
	//Total lattice/array dimensions.
	int xL = xU*nU;
	int yL = yU*mU;
	
	//Start lattice site (x,y); currently a cellular region.
	int xSP = (int)(xL/2.0) - (int)(xC/2.0) - 1;
	int ySP = (int)(yL/2.0) - (int)(yC/2.0) - 1;

	//'Modified' x,y-pos. of ith particle. Returns [0,xU-1 or yU-1]
	int modx;
	int mody;
	//Possible future position of particle. Returns: [0,xL-1 or yL-1]
	double testx = 0;
	double testy = 0;
	//'Modified' possible future position of particle.
	double modtestx = 0;
	double modtesty = 0;
	//State corresponding to particle position; (0: extracellular, 1: cellular).
	int state = 0;
	int teststate = 0;

	//Convert the ith particle position to an integer.
	int sxi;
	int syi;
	//Unit Cell Of Particle. x and y specify unique unit cell in 2D array.
	int ucopx;
	int ucopy;
	//Unit Cell Of Test [particle position].
	int ucotx;
	int ucoty;

	//Setting the start position of all particles. These arrays hold positions.
	double x[N] = {[0 ... (N-1)] = xSP};
	double y[N] = {[0 ... (N-1)] = ySP};

	//Particle density distribution 2D array; yL rows and xL columns.
	int rho[yL][xL];

	//Stepsize.
	double a = 1.0;
	//Probabilities for x or y motion.
	double px = 0.5, py = 0.5;
	//Stepping probabilities (arb. chosen) for intracellular.
	//Physical model: intracellular regions less diffusive (more viscous).
	double pnxi = 0.2, ppxi = 0.2;
	double pnyi = 0.2, ppyi = 0.2;
	//Stepping probabilities (arb. chosen) for extracellular.
	//Physical model: extracellular regions more diffusive (less viscous).
	double pnxe = 0.4, ppxe = 0.4;
	double pnye = 0.4, ppye = 0.4;
	//Coupled probabilities crossing from intra to extra (pie) or extra to intra (pei).
	//Although pie (or pei) arbitrarily chosen, these values are related. 
	double pie = 0.05;
	double pei = (pnxi/pnxe)*pie; 

	//Random variable and variables for MSD calculations.
	double rnd1, rnd2;
	double resultant_x, resultant_y, resultant_x_sq, resultant_y_sq, sum_resultant, msd;

	//File naming and preparation. Why do we get a segmentation fault below?
	//char *path = "/home/paul/Documents/thesis/particle-diffusion/data/";
	//char *f1 = strcat(path,"TEST1.txt");
	//char *f2 = strcat(path,"TEST1-stats.txt");
	char *f1 = "/home/paul/Documents/thesis/particle-diffusion/2D/2D-data/npr_test_3_xunit_500k.txt";
	char *f2 = "/home/paul/Documents/thesis/particle-diffusion/2D/2D-data/npr_test_3_xunit_500k_stats.txt";
	FILE *outdists, *outstats;
	outdists = fopen(f1, "w");
	outstats = fopen(f2, "w");
	
	//--------------------------------------------------------------------------
	//Below is where the main work is done.
	for (t = 0; t < tmax; t++){
		sum_resultant = 0;

		//Set density array elements to zero.
		for (j = 0; j < yL; j++){
			for(i = 0; i < xL; i++){
				rho[j][i] = 0;
			}
		}

		//Loop over all particles in system.
		for(n = 0; n < N; n++){

			//Determining the unit cell the particle is currently in.
			ucopx = (int)(x[n]/xU);
			ucopy = (int)(y[n]/yU);
			//Determining position of particle within unit cell.
			modx = x[n] - ucopx*xU;
			mody = y[n] - ucopy*yU;

			/*
			Usage of yC here does not reflect its meaning. Need to fix this.
			This is an issue if using more than one 'unit' cell in y-dir.
			I think the issue stems from having two different kinds of unit
			cells in the y-dir.
			*/
			if(mody < yC){
				/*
				Particle is in top half of 2D array and may be in either an
				intracellular or extracellular region.
				*/
				if(modx < xC){
					//Particle is currently in intracellular region.
					state = 1;
				}
				else{
					//Particle is currently in extracellular region.
					state = 0;
				}
			}
			else{
				/*
				The particle is in the lower half of the 2D array and is only
				in extracellular region. Refer to cell model.
				*/
				state = 0; //Maybe introduce a different state here?
			}
			if(state == 1){
				/*
				Apply intracellular conditions to test particle.
				Test membership of random number to 4 intervals.
				Each interval represents a different particle direction.
				*/
				rnd1 = (double)rand()/(double)RAND_MAX; //used for x motion.
				rnd2 = (double)rand()/(double)RAND_MAX; //used for y motion.

				/* 
				Updated algorithm; particle movement in x, followed by particle
				movement in y. Two moves per time step.
				*/

				//Motion in x.
				if(rnd1 < pnxi){
					//Generate position if to move -a in x-direction.
					testx = x[n] - a;
				}
				else if(rnd1 < pnxi + ppxi){
					//Generate position if to move +a in x-direction.
					testx = x[n] + a;
				}
				else{
					//Generate position if to stay in current position.
					testx = x[n];
				}

				//Motion in y.
				if(rnd2 < pnyi){
					//Generate position if to move -a in y-direction.
					testy = y[n] - a;
				}
				else if(rnd2 < pnyi + ppyi){
					//Generate position if to move +a in y-direction.
					testy = y[n] + a;
				}
				else{
					//Generate position if to stay in current position.
					testy = y[n];
				}

				/* THE FOLLOWING DOES NOT WORK AS INTENDED.
				Check event probability (events are independent).
				if(rnd1 < px){
					//Motion is in x direction. Position in y-direction fixed.
					testy = y[n];
					if(rnd2 < pnxi){
						//Generate position if to move -a in x-direction.
						testx = x[n] - a;
					}
					else if(rnd2 < pnxi + ppxi){
						//Generate position if to move +a in x-direction.
						testx = x[n] + a;
					}
					else{
						//Generate position if to stay in current position.
						testx = x[n];
					}
				}
				else{
					//Motion is in y direction. Position in x-direction fixed.
					testx = x[n];
					if(rnd2 < pnyi){
						//Generate position if to move -a in y-direction.
						testy = y[n] - a;
					}
					else if(rnd2 < pnyi + ppyi){
						//Generate position if to move +a in y-direction.
						testy = y[n] + a;
					}
					else{
						//Generate position if to stay in current position.
						testy = y[n];
					}
				}*/
			}
			else{
				/*
				Apply extracellular conditions to test particle.
				Test membership of random number to 4 intervals.
				Each interval represents a different particle direction.
				*/
				rnd1 = (double)rand()/(double)RAND_MAX; //used for x motion.
				rnd2 = (double)rand()/(double)RAND_MAX; //used for y motion.
				
				/* 
				Updated algorithm; particle movement in x, followed by particle
				movement in y. Two moves per time step.
				*/

				//Motion in x.
				if(rnd1 < pnxe){
					//Generate position if to move -a in x-direction.
					testx = x[n] - a;
				}
				else if(rnd1 < pnxe + ppxe){
					//Generate position if to move +a in x-direction.
					testx = x[n] + a;
				}
				else{
					//Generate position if to stay in current position.
					testx = x[n];
				}
				
				//Motion in y.
				if(rnd2 < pnye){
					//Generate position if to move -a in y-direction.
					testy = y[n] - a;
				}
				else if(rnd2 < pnye + ppye){
					//Generate position if to move +a in y-direction.
					testy = y[n] + a;
				}
				else{
					//Generate position if to stay in current position.
					testy = y[n];
				}
				
				/* THE FOLLOWING DOES NOT WORK AS INTENDED.
				Check event probability (events are independent).
				if(rnd1 < px){
					//Motion is in x direction. Position in y-direction fixed.
					testy = y[n];

					if(rnd2 < pnxe){
						//Generate position if to move -a in x-direction.
						testx = x[n] - a;
					}
					else if(rnd2 < pnxe + ppxe){
						//Generate position if to move +a in x-direction.
						testx = x[n] + a;
					}
					else{
						//Generate position if to stay in current position.
						testx = x[n];
					}
				}
				else{
					//Motion is in y direction. Position in x-direction fixed.
					testx = x[n];

					if(rnd2 < pnye){
						//Generate position if to move -a in y-direction.
						testy = y[n] - a;
					}
					else if(rnd2 < pnye + ppye){
						//Generate position if to move +a in y-direction.
						testy = y[n] + a;
					}
					else{
						//Generate position if to stay in current position.
						testy = y[n];
					}
				}*/
			}

			//Needed XOR logical operation.
			if(!(testx != x[n]) != !(testy != y[n])){
				/*
				Determine motion of the particle depending on test position.
				If particle is not set to move, nothing is done.
				*/
				if(testx >= 0 && testx < xL && testy >= 0 && testy < yL){
					/*
					If particle is within absolute array limits...
					Here is also determined if particle will cross into a new
					region. This is done by comparing the current and test state.
					*/
					ucotx = (int)(testx/xU);
					ucoty = (int)(testy/yU);
					modtestx = testx - ucotx*xU;
					modtesty = testy - ucoty*yU;

					//Usage of yC here does not reflect its meaning. Need to fix this.
					if(modtesty < yC){
						/*
						Particle is in top half of 2D array and may be in either an
						intracellular or extracellular region.
						*/
						if(modtestx < xC){
							//Particle is currently in intracellular region.
							teststate = 1;
						}
						else{
							//Particle is currently in extracellular region.
							teststate = 0;
						}
					}
					else{
						/*
						The particle is in the lower half of the 2D array and is only
						in extracellular region. Refer to cell model.
						*/
						teststate = 0; //Maybe introduce a different state here?
					}
					
					/*
					Has the particle moved to a different cell region?
					Particle positions are assigned in the next few if-else
					clauses depending on current state and future state. This is
					where particle positions are updated!
					*/
					if(state == teststate){
						/*
						Particle has not moved to different cell region. Simply
						set particle position to proposed position.
						*/
						x[n] = testx;
						y[n] = testy;
					}
					else if(state == 1 && teststate == 0){
						/*
						Particle is in intracellular region and at boundary.
						Future position is possibly extracellular region.
						Draw random number to determine movement.
						*/

						rnd1 = (double)rand()/(double)RAND_MAX;
						if(rnd1 < pie){
							//Particle will cross boundary into extracellular region.
							x[n] = testx;
							y[n] = testy;
						}
					}
					else{
						/*
						Particle is in extracellular region and at boundary.
						Future position is possibly intracellular region.
						Draw random number to determine movement.
						*/
						rnd1 = (double)rand()/(double)RAND_MAX;
						if(rnd1 < pei){
							//Particle will cross boundary in intracellular region.
							x[n] = testx;
							y[n] = testy;
						}					
					}
				}
			}
			//MSD calculations and file output below this line.-----------------
			resultant_x = x[n] - xSP;
			resultant_y = y[n] - ySP;
			resultant_x_sq = resultant_x*resultant_x;
			resultant_y_sq = resultant_y*resultant_y;

			sum_resultant += resultant_x_sq + resultant_y_sq;

			/*Convert particle position to integer and increment particle count 
			at that	position rho[sxi] by 1 unit.*/
			sxi = (int)(x[n]);
			syi = (int)(y[n]);
			rho[syi][sxi]++;
		}
		
		/*
		Writing density distribution data to file. Each time step is written as
		a new matrix with a space between each matrix. Each new line represents
		a row of data.
		*/
		for (j = 0; j < yL; j++){
			for(i = 0; i < xL; i++){
				fprintf(outdists, "%d ", rho[j][i]);
			}
			fprintf(outdists, "\n");
		}
		fprintf(outdists, "\n");

		//Writing mean-square-displacement data to file.
		//Done after looping over all particles; done once per time step.
		msd = sum_resultant/(double)N;
		fprintf(outstats, "%f\n", msd);
	}
	//return 0;
}

// get_argv(){
// 	//Get command-line arguments.
// 	//
// }