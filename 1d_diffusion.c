/* Diffusion simulation in 1D.*/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

/* CPP directives to define some constants used in the program.*/
#define N 10000 /*number of particles*/
#define L 1000  /*length of array*/
#define SP 500  /*start position of particles*/

int main(){
	int t, tmax = 5000;
	int i;

	float a = 1.0; /*stepsize*/
	float pl = 0.3, pr = 0.3, ps = 1.0-pl-pr; /*probablities of moving*/
	float rnd;
	float sum_x, sum_x2, avg_x, avg_x2;

	/* The following is an array to hold particle positions. 
	The ith element of that array cooresponds to the current 
	position of the ith particle. This information is never cleared. */
	float x[N] = {0};
	int sxi;

	int x_distribution[L] = {0};
	x_distribution[(int)SP-1] = N; /*setting the start pos. of particles*/

	char *f1 = "/home/paul/Documents/thesis/particle-diffusion/data/\
		1_diffusion_data_P030.txt";
	char *f2 = "/home/paul/Documents/thesis/particle-diffusion/data/\
		1d_diffusion_data_P030-stats.txt";
	FILE *dd, *outstats;
	dd = fopen(f1, "w");
	outstats = fopen(f2, "w");
	
	/*Writing in the initial ditribution...*/
	fprintf(dd, "#First integer is the time step.#\n");
		for(i = 0; i < L; i++){
		fprintf(dd, "%d ", x_distribution[i]);
	}
	fprintf(dd, "\n");


	for (t = 0; t < tmax; t++){
		/*For every time step...*/
		
		for (i = 0; i < L; i++){
			/*Loop over every cell in x_dist. array, set particle count to 0*/
			x_distribution[i] = 0;
		}
		sum_x = 0;
		sum_x2 = 0;
		
		for(i = 0; i < N; i++){
			/*Loop over all particles.*/
			sxi = (int)(x[i] + SP);
			
			/*Check the ith particle position and determine if at boundary.*/
			if(sxi != 0 && sxi != (L-1)){
				/*Particle is not at boundary, may move left, right, or stay.*/
				rnd = (float)rand()/(float)RAND_MAX;
				if(rnd < pl){
					/*The ith particle moves a distance 'a' to the left.*/
					x[i] -= a;
				}
				else if(rnd > pl && rnd < (pl + pr)){
					/*The ith particle moves a distance 'a' to the right.*/
					x[i] += a;
				}
			}
			else if(sxi == 0){
				rnd = (float)rand()/(float)RAND_MAX;
				if(rnd < 0.5){
					/*The ith particle moves a distance 'a' to the right.*/
					x[i] += a;
				}
			}
			else if(sxi == (L-1)){
				rnd = (float)rand()/(float)RAND_MAX;
				if(rnd < 0.5){
					/*The ith particle moves a distance 'a' to the left.*/
					x[i] -= a;
				}
			}
			/**/
			sum_x += x[i];
			sum_x2 += x[i]*x[i];

			sxi = (int)(x[i] + SP); /*Shift the position of the ith particle*/
			x_distribution[sxi]++; 
		}
		
		fprintf(dd, "%d ",t); /*first column is time*/
		for(i = 0; i < L; i++){
			fprintf(dd, "%d ", x_distribution[i]);
		}
		fprintf(dd,"\n");
		
		/*Writing some stats to file...*/
		avg_x = sum_x/(float)N;
		avg_x2 = sum_x2/(float)N;
		fprintf(outstats, "%f %f %f\n", avg_x, avg_x2, avg_x2-avg_x*avg_x);
	}
}