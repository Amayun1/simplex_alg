/*
Test program for the ez simplex method for linear programs in easy form.
Input is a linear program of the form
	mininmize c^Tx
	s.t.	  Ax <= b
		  x >=0
where b >= 0 and A is of full row-rank.

*/

// include simplex algorithm functions.
#include "simplex.h"
#include "simplex.c"

#include "print_matrix.c"

// Include time library for timing how long code takes to execute.
#include <time.h>

int main (void){

   // Test problem from Griva, Nash, Sofer (2009). See pages 97-98.
   lapack_int m = 3;
   lapack_int n = 2;

   double c[5] = {-1, -2};
   double A[6] = {-2, 1, -1, 1, 1, 0};
   double b[3] = {2, 3, 3};
   double x[2];

   // set up timing
   struct timespec start, end;

   timespec_get(&start, TIME_UTC); // start timer.

   double obj = ez_simplex(m, n, c, A, b, x);

   timespec_get(&end, TIME_UTC); // end timer.

   // calculate time in ms
   double time_spent = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_nsec - start.tv_nsec)/(1E+06);
   printf("EZ simplex finished in %f ms.\n", time_spent);

   printf("minimum x = \n");
   print_float_matrix(n, 1, x);

   printf("Objective value c^Tx = %f\n", obj);

   return 0;
}
