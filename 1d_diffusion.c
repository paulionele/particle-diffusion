/* Diffusion simulation in 1D.*/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

/* CPP directives to define some constants used in the program.*/
#define N 10000 /*number of particles*/
#define L 50  /*length of array*/
#define SP 25.0 /*start position of particles*/

int main(){
	int t, tmax = 500;
	int i;

	float a = 1.0; /*stepsize*/
	float pl = 0.3, pr = 0.3, ps = 1.0-pl-pr; /*probablities of moving*/
	float rnd;

	/*The following is an array to hold particle positions. 
	The ith element of that array cooresponds to the current 
	positon of the ith particle. This information is never cleared.*/
	float x[N] = {0};
	int sxi;

	int x_distribution[L] = {0};

	FILE *dd;
	dd = fopen("1d_diffusion_homogenous.txt", "w");

	for (t = 0; t < tmax; t++){
		/*For every time step...*/
		
		for (i = 0; i < L; i++){
			/*Loop over every cell in x_dist. array, set particle count to 0*/
			x_distribution[i] = 0;
		}

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

			sxi = (int)(x[i] + SP); /*Shift the position of the ith particle*/
			x_distribution[sxi]++; 
		}
		
		fprintf(dd, "%d",t);
		for(i = 0; i < L; i++){
			fprintf(dd, "%d ", x_distribution[i]);
		}
		fprintf(dd,"\n");
	}
}

/*
Some variables could be defined with the const keyword preceding, in the case
that we do not want to change the values of these variables later.
When you see const declaration in the function parameters, you know the 
function is not going to change the value of the parameter.
*/