/*

Function for performing the simplex algorithm in C.
Requires linking lpacke and cblas to use:

-llpacke -lapack -lblas.


Author: Charles Panigeo
Date Created: 6/24/2026
*/

#include "simplex.h"

void sf_simplex_phase_one(lapack_int m, lapack_int n, const double *A, const double *rhs, volatile double *bfs){
/*
Calculcates an initial BFS of the linear program in standard form Ax = b using the phase-one method.

Inputs:
   lapack_int m: Number of rows of the contraints matrix A and rhs b.
   lapack_int n: Number of columns of A and the number of variables to optimize.
   const double *A: Pointer to an array of length m * n for the constaint matrix. Must be in row major order.
   const double *rhs: Pointer to an array of length m for the rhs of the constraints.
   volatile double *bfs: Pointer to an array of length n where the calculated BFS should be stored.
Outputs:
   x (overwrites bfs): Coordinates for a bfs of Ax = b.
*/

   // Create modifed constraints matrix.
   // From A, create size for m * (n + m) matrix [A I_m].
   double *temp_A = malloc(m * (n + m) * sizeof(double));
   if(temp_A == NULL){
      printf("Memory allocation failed! Aborting function %s\n", __func__);
      exit(EXIT_FAILURE);
   }

   // Copy in A to the left n columns of of [A I_m]
   for(int rows = 0; rows < m; rows++){
      for(int cols = 0; cols < n; cols++){
         temp_A[(n + m) * rows + cols] = A[n * rows + cols];
      }
   }
   // Set right m columns to be I_m
   for(int rows = 0; rows < m; rows++){
      for(int cols = 0; cols < m; cols++){
         if(cols == rows){
            temp_A[n + (n + m) * rows + cols] = 1;
         }
         else{
            temp_A[n + (n + m) * rows + cols] = 0;
         }
      }
   }

   // create c^T for the phase-one problem. Array of n + m with zeros in first n columns and 1 for the artificial variables.
   double *c = calloc((n + m), sizeof(double));
   if(c == NULL){
      printf("Memory allocation failed! Aborting function %s\n", __func__);
      exit(EXIT_FAILURE);
   }
   for(int i = n; i < n + m; i++){
      c[i] = 1;
   }

   // create BFS for the phase-one problem. 0's in first n columns, 1s for each artificial variable
   double *x = calloc((n + m), sizeof(double));
   if(x == NULL){
      printf("Memory allocation failed! Aborting function %s\n", __func__);
      exit(EXIT_FAILURE);
   }
   for(int i = 0; i < m; i++){
      x[i + n] = rhs[i];
   }

    // Run simplex with Ax = b and c^Tx.
   double object = sf_simplex(m, m + n, c, temp_A, rhs, x);

   // Check if objective value is zero to make sure the artificial variables are 0;
   if(object != 0){
      printf("Phase-one simplex failed! Calculated BFS has artificial variables. Aborting program.");
      exit(EXIT_FAILURE);
   }

   // Write first n coordinates to bfs, then double check it is a bfs.
   for(int i = 0; i < n; i++){
      bfs[i] = x[i];
   }
   if(!is_basic(m, n, A, rhs, bfs)){
      printf("Phase-one simplex failed! Calculated solution is not a BFS of original proble. Aborting program.");
      exit(EXIT_FAILURE);
   }

   // Free allocated memory and clean dangling pointers.
   free(temp_A);
   temp_A = NULL;
   free(c);
   c == NULL;
   free(x);
   x == NULL;

   return;
}

double sf_simplex(lapack_int m, lapack_int n, const double *c, const double *A, const double *rhs, volatile double *bfs){
/*
Program which implements the simplex method for linear programs in standard form.
Input is a linear program of the form
	mininmize c^Tx
	s.t.	  Ax = b
		  x >=0
where b >= 0 and A is of full row-rank.

Inputs:
   lapack_int m: Number of rows of the contraints matrix A and rhs b.
   lapack_int n: Number of columns of A and the number of variables to optimize.
   const double *c: Pointer to an array length n containing the weights of the objective function.
   const double *A: Pointer to an array of length m*n for the constaint matrix. Must be in row major order.
   const double *rhs: Pointer to an array of length m for the rhs of the constraints.
   volatile double *bfs: Pointer to an array of length n containing the starting BFS.
Outputs:
   x (overwrites bfs): Coordinates for an optimal solution
   return value: Value of c^Tx for the calculated minimum.
*/

   // checks that b >= 0.
   for(int i = 0; i < m; i++){
      if(rhs[i] < 0){
         printf("Problem is not in standard form. Aborting program. Check rhs of Ax=b and try again.\n");
         exit(EXIT_FAILURE);
      }
   }

   // check if feasible.
   for(int i = 0; i < n; i++){
      if(bfs[i] < 0){
         printf("Warning: Given BFS is not feasible. Proceed anyway? Enter Y to override warning. \n");

         // Gets one char from user. It it is Y or y, proceed anyway. If anything else, abort program.
         char ch;
         scanf(" %c", &ch);
         if(tolower(ch) == 'y'){
            printf("Overriding warning! \n");
            break;
         }
         else{
            printf("Aborting program. Check BFS before trying again.\n");
            exit(EXIT_FAILURE);
         }
      }
   }

   // Check if close to a basic solution.
   //bool is_basic(int m, int n, const double *matrix, const double *rhs, const double *bfs)

   if(is_basic(m, n, A, rhs, bfs) != true){
      printf("Warning: Given BFS is not close to a basic soluion. Proceed anyway? Enter Y to override warning. \n");

      // Gets one char from user. It it is Y or y, proceed anyway. If anything else, abort program.
      char ch;
      scanf(" %c", &ch);
      if(tolower(ch) == 'y'){
         printf("Overriding warning! \n");
      }
      else{
         printf("Aborting program. Check BFS before trying again.\n");
         exit(EXIT_FAILURE);
      }
   }

   // From starting BFS, initalize BFS indices array.
   bool *BFS_indices = malloc(n * sizeof(bool));

   if (BFS_indices == NULL){
      printf("Memory allocation failed! Aborting function %s\n", __func__);
      exit(EXIT_FAILURE);
   }

   for(int i = 0; i < n; i++){
      if(bfs[i] > 0){
         BFS_indices[i] = true;
      }
      else{
         BFS_indices[i] = false;
      }
   }

   // allocate memory for B and N matrices, c_B, c_N, A_t/hatA_t, b/hatb.
   double *B = malloc(m * m * sizeof(double));
   if (B == NULL){
      printf("Memory allocation failed! Aborting function %s\n", __func__);
      exit(EXIT_FAILURE);
   }

   double *N = malloc(m * (n - m) * sizeof(double));
   if (N == NULL){
      printf("Memory allocation failed! Aborting function %s\n", __func__);
      exit(EXIT_FAILURE);
   }

   double *c_B = malloc(m * sizeof(double));
   if (c_B == NULL){
      printf("Memory allocation failed! Aborting function %s\n", __func__);
      exit(EXIT_FAILURE);
   }

   double *c_N = malloc((n - m) * sizeof(double));
   if (c_N == NULL){
      printf("Memory allocation failed! Aborting function %s\n", __func__);
      exit(EXIT_FAILURE);
   }

   double *A_t = malloc(m * sizeof(double));
   if (A_t == NULL){
      printf("Memory allocation failed! Aborting function %s\n", __func__);
      exit(EXIT_FAILURE);
   }

   double *b = malloc(m * sizeof(double));
   if (b == NULL){
      printf("Memory allocation failed! Aborting function %s\n", __func__);
      exit(EXIT_FAILURE);
   }

   // Allocate memory for ipiv and info for LU factorization.
   lapack_int *ipiv = (lapack_int*) malloc(m * m * sizeof(lapack_int));
   if(ipiv == NULL){
      printf("Memory allocation failed! Aborting function %s\n", __func__);
      exit(EXIT_FAILURE);
   }

   // MAIN SIMPLEX LOOP
   start_simplex:
   int i = 0;
   while(i < MAX_SIMPLEX_ITERS){

      // Copy in B matrix from BFS indices. cursor_B and cursor_N are "cursor variables" for the columns of B and N respectively.
      // The loop runs over each column of A and decides based on the BFS indices whether or not to copy the values
      // into the next column of either B or N. Cursors are initialized at 0.
      // Also copies the corrsponding objective values from c into c_B and c_N.
      int cursor_B = 0, cursor_N = 0;
      for(int j = 0; j < n; j++){
         if(BFS_indices[j] == true){
            // If the current BFS index is true, copy current column into next column of B and increment the B cursor.
            for(int i = 0; i < m; i++){
               B[m * i + cursor_B] = A[n * i + j];
            }
            c_B[cursor_B] = c[j];
            cursor_B++;
         }
         else{
            // Else, copy the column into the next column of N instead and increment the N cursor.
            for(int i = 0; i < m; i++){
               N[(n - m) * i + cursor_N] = A[n * i + j];
            }
            c_N[cursor_N] = c[j];
            cursor_N++;
         }
      }

      // Calculate LU factorization for B. The result overwrites B.
      lapack_int info;
      info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, m, B, m, ipiv);

      if (info != 0){
         printf("LU Factorization failed! Aborting function %s.\n", __func__);
         exit(EXIT_FAILURE);
      }

      // Solve Bx_b = b, overwriting b with x_b = hatb
      // copy rhs values into b for solving.
      memcpy(b, rhs, m * sizeof(double));
      info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', m, 1, B, m, ipiv, b, 1);

      if(info != 0){
         printf("Error! Linear system not solvable in call of %s.\n", __func__);
         exit(EXIT_FAILURE);
      }

      // calculate simplex multipliers B^Ty = c_B, overwriting y into c_B
      info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'T', m, 1, B, m, ipiv, c_B, 1);
      if(info != 0){
         printf("Error! Linear system not solvable in call of %s.\n", __func__);
         exit(EXIT_FAILURE);
      }

      // Calculate hatc_n = - N^Ty + c_N, overwriting hatc_n into c_N.
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n - m, 1, m, -1, N, n - m, c_B, 1, 1, c_N, 1);

      // Test for optimality. Finds min value of c_N, noting it's index in c_N and check if >= 0)
      double min = INFINITY;
      int index_min = 0;
      for(int i = 0; i < n - m; i++){
         if(c_N[i] < min){
            min = c_N[i];
            index_min = i;
         }
      }

      // If optimal, break the simplex loop.
      if(min >= 0){
         break;
      }

      // Set entering variable for the next loop.
      // Based on the index_min from c_N, find the correponding index in BFS_indices.
      // From left to right, we keep track of how many non_BFS indices we've seen so far indcrementing the seen variable
      // Then checking if we've seen index_min non BFS variables yet of them
      // Once we get to the index_min'th one, we change it to true as the entering variable.
      int seen = 0;
      int t; //entering variable
      for(int i = 0; i < n; i++){
         if(BFS_indices[i] == false){
            if(seen == index_min){
               t = i;
               break;
            }
            seen++;
         }
      }

      // Calculate the leaving variable.

      // copy column t of A into A_t
      for(int i = 0; i < m; i++){
         A_t[i] = A[n * i + t];
      }

      // Solve B hatA_t = A_t, overwriting solution into A_t.
      info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', m, 1, B, m, ipiv, A_t, 1);
      if(info != 0){
         printf("Error! Linear system not solvable in call of %s.\n", __func__);
         exit(EXIT_FAILURE);
      }

      // If hata_t > 0, calculate hatb/hata_t, store, and increment num_pos. Reusing c_B storing.
      // Else store infinity instead.
      int num_pos = 0;
      for(int i = 0; i < m; i++){
         if(A_t[i] > 0){
            c_B[i] = b[i]/A_t[i];
            num_pos++;
         }
         else{
            c_B[i] = INFINITY;
         }
      }

      // Check if problem is unbounded.
      if(num_pos == 0){
         printf("Linear program is unbounded! Aborting. \n");
         exit(EXIT_FAILURE);
      }

      // calculate leaving variable index A_t (index among non-BFS variables).
      // reusing min and index_min from earlier.
      min = INFINITY;
      index_min = 0;
      for(int i = 0; i < m; i++){
         if(c_B[i] < min){
            min = c_B[i];
            index_min = i;
         }
      }

      // Set leaving variable for the next loop.
      // Based on the index_min, find the correponding index in BFS_indices.
      // From left to right, we keep track of how many BFS indices we've seen so far incrementing the seen variable
      // Then checking if we've seen index_min number of BFS variables yet of them
      // Once we get to the index_min'th one, we change it to false as the leaving variable.
      seen = 0;
      int s; //leaving variable
      for(int i = 0; i < n; i++){
         if(BFS_indices[i]){
            if(seen == index_min){
               s = i;
               break;
            }
            seen++;
         }
      }

      // Update entering and leaving variables
      BFS_indices[t] = true;
      BFS_indices[s] = false;

      i++; //increment number of simplex steps
   }

   // Take hat and convert to coordinates, overwriting BFS.
   // Cursor points to the next hatb which is yet to be copied.
   int cursor = 0;
   for(int i = 0; i < n; i++){
      // if the index is active, copy in the coresp entry from hat and advance cursor.
      if(BFS_indices[i] == true){
         bfs[i] = b[cursor];
         cursor++;
      }
      // if the index is not active, write 0 and don't advance.
      else{
         bfs[i] = 0;
      }
   }

   // Tests if MAX_ITERS was reached. If so, the value in bfs is not optimal.
   // Prompt user if they want to run another round of MAX_ITERS iterations.
   if(i == MAX_SIMPLEX_ITERS){
      printf("Simplex reached the max allowed iterations without finding an optimal solution. \n");
      printf("Try again with another %d iterations? (Y/n):\n", MAX_SIMPLEX_ITERS);

      // Gets one char from user. It it is Y or y, proceed anyway. If anything else, abort program.
      char ch;
      scanf(" %c", &ch);
      if(tolower(ch) == 'y'){
         printf("Running simplex again! \n");
         goto start_simplex;
      }
      else{
         printf("Ending simplex without finding the optimal solution. Proceed with caution. \n");
      }
   }

   // Free allocated memory blocks and clean up dangling pointers.
   free(BFS_indices);
   BFS_indices = NULL;
   free(B);
   B = NULL;
   free(N);
   N = NULL;
   free(c_B);
   c_B = NULL;
   free(c_N);
   c_N = NULL;
   free(A_t);
   A_t = NULL;
   free(b);
   b = NULL;
   free(ipiv);
   ipiv = NULL;

   // calculate optimal value and return.
   return cblas_ddot(n, bfs, 1, c, 1);
}

double ez_simplex(lapack_int m, lapack_int n, const double *c, const double *A, const double *rhs, volatile double *x){
/*
Wrapper for sf_simplex.
Solves a linear program of the form
	mininmize c^Tx
	s.t.	  Ax <= b
		  x >=0
where b >= 0 and A is of full row-rank.

The problem is converted to standard form by adding slack variables to every constraint.
The initial BFS uses every slack variable, taking [0, b].

Inputs:
   lapack_int m: Number of rows of the contraints matrix A and rhs b.
   lapack_int n: Number of columns of A and the number of variables to optimize.
   const double *c: Pointer to an array length n containing the weights of the objective function.
   const double *A: Pointer to an array of length m*n for the constaint matrix. Must be in row major order.
   const double *rhs: Pointer to an array of length m for the rhs of the constraints.
   volatile double *x: Pointer to an array of length n where the optimal solution is stored.
Outputs:
   x: Coordinates for an optimal solution
   return value: Value of c^Tx for the calculated minimum.
*/

   // Create modifed constraints matrix.
   // From A, create size for m * (n + m) matrix [A I_m].
   double *temp_A = malloc(m * (n + m) * sizeof(double));
   if(temp_A == NULL){
      printf("Memory allocation failed! Aborting function %s\n", __func__);
      exit(EXIT_FAILURE);
   }

   // Copy in A to the left n columns of of [A I_m]
   for(int rows = 0; rows < m; rows++){
      for(int cols = 0; cols < n; cols++){
         temp_A[(n + m) * rows + cols] = A[n * rows + cols];
      }
   }
   // Set right m columns to be I_m
   for(int rows = 0; rows < m; rows++){
      for(int cols = 0; cols < m; cols++){
         if(cols == rows){
            temp_A[n + (n + m) * rows + cols] = 1;
         }
         else{
            temp_A[n + (n + m) * rows + cols] = 0;
         }
      }
   }

   // create c^T with m more zeroes as the new objective function.
   double *temp_c = calloc(n + m, sizeof(double));
   if(temp_c == NULL){
      printf("Memory allocation failed! Aborting function %s\n", __func__);
      exit(EXIT_FAILURE);
   }
   for(int i = 0; i < n; i++){
      temp_c[i] = c[i];
   }

   // create BFS for the problem. 0's in first n columns, rhs for each slack variable position.
   double *bfs = calloc((n + m), sizeof(double));
   if(bfs == NULL){
      printf("Memory allocation failed! Aborting function %s\n", __func__);
      exit(EXIT_FAILURE);
   }
   for(int i = 0; i < m; i++){
      bfs[n + i] = rhs[i];
   }

   // pass problem into standard form simplex alg.
   double objective = sf_simplex(m, n + m, temp_c, temp_A, rhs, bfs);

   // copy solution into solution pointer without slack variables.
   for(int i = 0; i < n; i++){
      x[i] = bfs[i];
   }

   // Free allocated memory blocks and clean dangling pointers.
   free(temp_A);
   temp_A = NULL;
   free(temp_c);
   temp_c = NULL;
   free(bfs);
   bfs = NULL;

   // Return objective value.
   return objective;
}

bool is_basic(int m, int n, const double *matrix, const double *rhs, const double *bfs){
   // Checks if bfs really is a basic solution to Ax = b, where A is an m * n matrix.
   // Returns true or false.

   // create temp variable the calculated Ax - b.
   double *solution = malloc(n * sizeof(double));
   if(solution == NULL){
      printf("Memory allocation failed! Aborting function %s.\n", __func__);
      exit(EXIT_FAILURE);
   }

   // Copy rhs values into solution
   memcpy(solution, rhs, m * sizeof(double));

   // Calculate Ax - b and overwrite b.
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, 1, n, 1.0, matrix, n, bfs, 1, -1.0, solution, 1);

   bool is_basic_flag;
   // check if the 2-norm of solution is large compared to the 2-norm of rhs.
   if(cblas_dnrm2(m, solution, 1) > 1.0E-06 * cblas_dnrm2(m, rhs, 1)){
      is_basic_flag = false;
   }
   else{
      is_basic_flag = true;
   }

   // free allocated memory and return.
   free(solution);
   solution = NULL;

   return is_basic_flag;
}
