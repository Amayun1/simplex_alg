/*
Test program for the ez simplex method for linear programs in easy form.
Input is a linear program of the form
	mininmize c^Tx
	s.t.	  Ax <= b
		  x >=0
where b >= 0 and A is of full row-rank.

The test problem is the Klee-Minty cube.
https://en.wikipedia.org/wiki/Klee%E2%80%93Minty_cube

*/

#define DIMENSIONS 26 // Number of dimensions for the cube.

//Override max simplex iters.
#define MAX_SIMPLEX_ITERS (1UL << (DIMENSIONS))


// include simplex algorithm functions.
#include "simplex.h"
#include "simplex.c"

// Include time library for timing how long code takes to execute.
#include <time.h>

#include "print_matrix.c"

int main (void){

   printf("Setting up the linear program!\n");

   lapack_int m = DIMENSIONS;
   lapack_int n = DIMENSIONS;

   // create a look-up table for pows of 2 and 5.
   // math.h doesn't work since it uses floating-point and Klee-Minty is very sensitive to the constraints.

   long long pows_of_2[DIMENSIONS + 1] = {1};
   long long pows_of_5[DIMENSIONS + 1] = {1};
   for(int i = 1; i <= DIMENSIONS; i++){
      pows_of_2[i] = 2 * pows_of_2[i - 1];
      pows_of_5[i] = 5 * pows_of_5[i - 1];
   }

   // Create the objective function
   double c[n];
   for(int i = 1; i <= n; i++){
      c[i] = - (double) pows_of_2[n - i];
   }
   // Create rhs of contstraints
   double b[n];
   for(int i = 0; i < n; i++){
      b[i] = (double) pows_of_5[i + 1];
   }

   // create contraints matrix
   double A[n * n];
   for(int rows = 0; rows < n; rows++){
      for(int cols = 0; cols < n; cols++){
         if(cols < rows){
            A[n * rows + cols] = (double) pows_of_2[rows - cols + 1];
         }
         else{
            A[n * rows + cols] = 0;
         }
      }
      A[n * rows + rows] = 1;
   }

   // Set up array to store the solution.
   double x[n];

   printf("Set up complete. Starting the calculation.\n");

   // set up timing
   struct timespec start, end;

   timespec_get(&start, TIME_UTC); // start timer.

   double obj = ez_simplex(m, n, c, A, b, x);

   timespec_get(&end, TIME_UTC); // end timer.

   // calculate time in ms
   double time_spent = (end.tv_sec - start.tv_sec) * + (end.tv_nsec - start.tv_nsec)/(1E+09);
   printf("EZ simplex finished in %f s.\n", time_spent);

   printf("minimum x = \n");
   print_float_matrix(n, 1, x);

   printf("Objective value c^Tx = %f\n", obj);

   return 0;
}
