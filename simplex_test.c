/*
Program which implements the simplex method for linear programs in standard form.
Input is a linear program of the form
	mininmize c^Tx
	s.t.	  Ax = b
		  x >=0
where b >= 0 and A is of full row-rank.

*/

// include simplex algorithm functions.
#include "simplex.h"
#include "simplex.c"

// Include time library for timing how long code takes to execute.
#include <time.h>

int main (void){

   // set up timing
   struct timespec start, end;

   // Test problem from Griva, Nash, Sofer (2009). See pages 132-133.
   lapack_int m = 3;
   lapack_int n = 5;

   double c[5] = {-1, -2, 0, 0, 0};
   double A[15] = {-2, 1, 1, 0, 0, -1, 2, 0, 1, 0, 1, 0, 0, 0, 1};
   double b[3] = {2, 7, 3};
   double bfs[5];

   // call simplex_phase_one to get a starting bfs for simplex method.
   sf_simplex_phase_one(m, n, A, b, bfs);
   printf("Calculated BFS = \n");
   print_float_matrix(n, 1, bfs);

   timespec_get(&start, TIME_UTC); // start timer.

   double obj = sf_simplex(m, n, c, A, b, bfs);

   timespec_get(&end, TIME_UTC); // end timer.

   // calculate time in ms
   double time_spent = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_nsec - start.tv_nsec)/(1E+06);
   printf("Simplex finished in %f ms.\n", time_spent);

   printf("minimum x = \n");
   print_float_matrix(n, 1, bfs);

   printf("Objective value c^Tx = %f\n", obj);


   

   return 0;
}
