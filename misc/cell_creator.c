/* 
This script is not functioning as expected. An unresolved issue of pre-determining
the length of the total inhomogenous array exists. The length of the array needs
to be known and defined through a CPP directive, otherwise a compile-time error
results. Maybe there is a workaround? AFAIK, everything else works.

All integer numbers used are for demonstration purposes and printf statements
are for poor mans debugging.
*/

#include "stdio.h"
#include "stdlib.h"

/*Total array length including all subunits.
  Need to this define this up here who the f knows.*/
#define n 100???

int main(){
	/*For loop variables counter variables.*/
	int a;
	int b;
	int c;

	/*Counter to deal with problems of mid-level language lol.*/
	int d = 0;

	/*Flip the j,k bits for change starting with intra or extracellular.*/
	int j = 1;
	int k = 0;

	int lic = 4;/*number of cells in intracellular super cell*/
	int lec = 2;/*number of cells in extracellular super cell*/
	int lsc = 7;/*number of super cells, this is okay as the for-loop limit*/

	/*Generated array, the total length including all subcells
	  needs to be known first and likely initialized in the CPP define directions.*/
	int ct[n] = {0};

	for(a = 0; a < lsc; a++){
		if(j!=0){
			for(b = 0; b < lic; b++){
				d++;
				ct[d] = 1;
			}
			j = 0;
			k = 1;
			printf("%s\n","hit1");
		}
		else if (k!=0){
			for(c = 0; c < lec; c++){
				d++;
				ct[d] = 0;
			}
			j = 1;
			k = 0;
			printf("%s\n","hit2");
		}
		printf("%s\n","hello");
	}
	int z;
	for(z = 0; z < lsc; z++){printf("%d ",ct[z]);}
	printf("\n");
}