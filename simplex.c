/*

Function for performing the simplex algorithm in C.
Requires linking LAPACKE, CBLAS, and math to use:

-llpacke -lapack -lblas -lm.


Author: Charles Panigeo
Date Created: 6/24/2026
*/

#include "simplex.h"

void sf_simplex_phase_one(lapack_int m, lapack_int n, const double *A, const double *rhs, volatile double *bfs, volatile bool *BFS_indices){
/*
Calculcates an initial BFS of the linear program in standard form Ax = b using the phase-one method.

Inputs:
   lapack_int m: Number of rows of the contraints matrix A and rhs b.
   lapack_int n: Number of columns of A and the number of variables to optimize.
   const double *A: Pointer to an array of length m * n for the constaint matrix. Must be in row major order.
   const double *rhs: Pointer to an array of length m for the rhs of the constraints.
   volatile double *bfs: Pointer to an array of length n where the calculated BFS should be stored.
   volatile bool *BFS_indices: Array where the BFS indices for the calculates solution should be stored.
Outputs:
   x (overwrites bfs): Coordinates for a bfs of Ax = b.
   BFS_indices: BFS indices for the calculated iniital BFS.
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

   // create temp BFS indicies for the phase-one problem.
   bool *temp_BFS_indices = malloc((n + m) * sizeof(bool));
      if(temp_BFS_indices == NULL){
      printf("Memory allocation failed! Aborting function %s\n", __func__);
      exit(EXIT_FAILURE);
   }

   //initialize temp_BFS_indices with the BFS for the phase-one problem.
   for(int i = 0; i < m + n; i++){
      if(i < n){
         temp_BFS_indices[i] = false;
      }
      else{
         temp_BFS_indices[i] = true;
      }
   }

   // Run simplex with Ax = b and c^Tx.
   double object = sf_simplex(m, m + n, c, temp_A, rhs, x, temp_BFS_indices);

   // Check if objective value is zero to make sure the artificial variables are 0;
   if(object > 0){
      printf("Phase-one simplex failed! Calculated BFS has artificial variables. Aborting program. \n");
      exit(EXIT_FAILURE);
   }

   // Write first n coordinates to bfs, then double check it is a bfs.
   for(int i = 0; i < n; i++){
      bfs[i] = x[i];
      BFS_indices[i] = temp_BFS_indices[i];
   }
   if(!is_basic(m, n, A, rhs, bfs)){
      printf("Phase-one simplex failed! Calculated solution is not a BFS of original problem. Aborting program.");
      exit(EXIT_FAILURE);
   }

   // Free allocated memory and clean dangling pointers.
   free(temp_A);
   temp_A = NULL;
   free(c);
   c = NULL;
   free(x);
   x = NULL;
   free(temp_BFS_indices);
   temp_BFS_indices = NULL;

   return;
}

double sf_simplex(lapack_int m, lapack_int n, const double *c, const double *A, const double *rhs, volatile double *bfs, volatile bool *BFS_indices){
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
   volatile bool *BFS_indices: pointer to an array of bools which specify which indices are active for the starting BFS.
Outputs:
   x (overwrites bfs): Coordinates for an optimal solution
   return value: Value of c^Tx for the calculated minimum.
   BFS_indices: Contains the bfs indices for the optimal solution.
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

   // info variable for LU factorization.
   lapack_int info;

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

   // create BFS_indices for the supplied BFS.
   bool *BFS_indices = malloc((n + m) * sizeof(bool));
   if(BFS_indices == NULL){
      printf("Memory allocation failed! Aborting function %s\n", __func__);
      exit(EXIT_FAILURE);
   }
   for(int i = 0; i < m + n; i++){
      if(i < n){
         BFS_indices[i] = false;
      }
      else{
         BFS_indices[i] = true;
      }
   }

   // pass problem into standard form simplex alg.
   double objective = sf_simplex(m, n + m, temp_c, temp_A, rhs, bfs, BFS_indices);

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
   free(BFS_indices);
   BFS_indices = NULL;

   // Return objective value.
   return objective;
}

bool is_basic(int m, int n, const double *A, const double *rhs, const double *bfs){
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
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, 1, n, 1.0, A, n, bfs, 1, -1.0, solution, 1);

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

bool sf_verify_sol(lapack_int m, lapack_int n, const double *c, const double *A, const double *rhs, const double *x){
/*
Verifies if a given vector x is a solution to a linear program in standard form
   min c^Tx
   s.t. Ax = b
   x >=0, b >=0.
The function solves the dual linear program
   max b^Ty
   s.t. A^Ty <= c
and checks if c^Tx = b^Ty. The dual LP is converted to standard form to call sf_simplex.

Inputs:
   lapack_int m: Number of rows of the contraints matrix A and rhs b.
   lapack_int n: Number of columns of A and the number of variables to optimize.
   const double *c: Pointer to an array length n containing the weights of the objective function.
   const double *A: Pointer to an array of length m*n for the constaint matrix. Must be in row major order.
   const double *rhs: Pointer to an array of length m for the rhs of the constraints.
   const double *x: Pointer to an array of length n where the proposed optimal solution is stored.
Outputs:
   bool return value: True if x is optimal, false if x is not optimal.
*/

   // Allocate memory for the dual LP problem parameters (converted into standard form).
   double *temp_A_trans = malloc(n * (m + m + n) * sizeof(double));
   if(temp_A_trans == NULL){
      printf("Memory allocation failed! Aborting function %s.\n", __func__);
      exit(EXIT_FAILURE);
   }

   double *temp_c = malloc(n * sizeof(double));
   if(temp_c == NULL){
      printf("Memory allocation failed! Aborting function %s.\n", __func__);
      exit(EXIT_FAILURE);
   }

   // construct temp_c, taking absolute values so the RHS of temp_A^Ty = c is nonnegative.
   // also add more zeroes so it is of the full m + m + n length.
   for(int i = 0; i < n; i++){
      temp_c[i] = fabs(c[i]);
   }

   // construct contraints matrix for the dual LP in standard form.
   // we replace each variable y_i = y_i' - y_i'' so we can have nonnegativity contraints on each variable.
   // [A^T -A^T I] with jth row multiplied by sign(c_j) as necessary.
   for(int rows = 0; rows < n; rows++){
      for(int cols = 0; cols < m; cols++){
         //need to check if c[rows] is negative to see if we need to multiply the row by -1.
         if(c[rows] >= 0){
            temp_A_trans[rows * (m + m + n) + cols] = A[n * cols + rows]; // copy in A^T
            temp_A_trans[rows * (m + m + n) + m + cols] = -A[n * cols + rows]; // copy in -A^T
         }
         else{
            temp_A_trans[rows * (m + m + n) + cols] = -A[n * cols + rows]; // copy in -A^T
            temp_A_trans[rows * (m + m + n) + m + cols] = A[n * cols + rows]; // copy in A^T
         }
      }
      for(int cols = 0; cols < n; cols++){
         // copy in -I (with row multiplied by -1 if -c[rows] is negative if necessary)
         if(rows == cols){
            if(c[rows] >=0){
               temp_A_trans[rows * (m + m + n) + m + m + cols] = 1;
            }
            else{
               temp_A_trans[rows * (m + m + n) + m + m + cols] = -1;
            }
         }
         else{
            temp_A_trans[rows * (m + m + n) + m + m + cols] = 0;
         }
      }
   }

   // allocate memory for the BFS/solution of the modified dual LP
   double *y = malloc((m + m + n) * sizeof(double));
   if(y == NULL){
      printf("Memory allocation failed! Aborting function %s.\n", __func__);
      exit(EXIT_FAILURE);
   }

   // Allocate memory for the BFS_indices of the modified dual LP
   bool *BFS_indices = malloc((m + m + n) * sizeof(bool));
   if(BFS_indices == NULL){
      printf("Memory allocation failed! Aborting function %s.\n", __func__);
      exit(EXIT_FAILURE);
   }

   // call sf_simplex_phase_one to get a BFS
   sf_simplex_phase_one(n, m + m + n, temp_A_trans, temp_c, y, BFS_indices);

   // allocate memory for new objective function since there are now m + n more variables in it.
   // also multiply by -1 to convert max into min for standard form.
   double *new_obj = calloc(m + m + n, sizeof(double));
   if(new_obj == NULL){
      printf("Memory allocation failed! Aborting function %s.\n", __func__);
      exit(EXIT_FAILURE);
   }

   // new objective function
   for(int i = 0; i < m; i++){
      new_obj[i] = -rhs[i];
      new_obj[i + m] = rhs[i];
   }

   // call sf_simplex to solve the modified dual LP
   sf_simplex(n, m + m + n, new_obj, temp_A_trans, temp_c, y, BFS_indices);

   // convert the modified y back into the solution for the original dual LP.
   // stores this back in to the same array so we don't need to allocate another array.
   // The entries beyond the first n are now garbage, but we never access them again so this is okay.
   for(int i = 0; i < m; i++){
      y[i] = y[i] - y[i + m];
   }

   printf("Primal-dual gap = %f.\n", fabs(cblas_ddot(n, c, 1, x, 1) - cblas_ddot(m, rhs, 1, y, 1)));

   bool is_solution_flag;
   // test if c^Tx - b^Ty > tolerance.
   // The choice of tolerance is equite arbitrary, but if c and b are huge, I think this method is a more "fair" comparison.
   if(fabs(cblas_ddot(n, c, 1, x, 1) - cblas_ddot(m, rhs, 1, y, 1))
       < 1E-06 * (cblas_dnrm2(n, c, 1) + cblas_dnrm2(m, rhs, 1)) / 2){
      is_solution_flag = true;
   }
   else{
      is_solution_flag = false;
   }

   // Free allocated memory and clean dangling pointers
   free(temp_A_trans);
   temp_A_trans = NULL;
   free(temp_c);
   temp_c = NULL;
   free(y);
   y = NULL;
   free(BFS_indices);
   BFS_indices = NULL;

   return is_solution_flag;
}

// TODO write function to verify a linear program in EZ form by converting to standard form and invoking sf_verify_sol.

bool ez_verify_sol(lapack_int m, lapack_int n, const double *c, const double *A, const double *rhs, const double *x){
/*
Verifies if a given vector x is a solution to a linear program in EZ form
   min c^Tx
   s.t. Ax <= b
   x >=0, b >=0.
Converts the problem and the proposed optimal solution to standard form and calls sf_verify_sol.

Inputs:
   lapack_int m: Number of rows of the contraints matrix A and rhs b.
   lapack_int n: Number of columns of A and the number of variables to optimize.
   const double *c: Pointer to an array length n containing the weights of the objective function.
   const double *A: Pointer to an array of length m*n for the constaint matrix. Must be in row major order.
   const double *rhs: Pointer to an array of length m for the rhs of the constraints.
   const double *x: Pointer to an array of length n where the proposed optimal solution is stored.
Outputs:
   bool return value: True if x is optimal, false if x is not optimal.
*/

   // Allocate memory for standard form matrix [A I] = temp_A
   double *temp_A = malloc(m * (m + n) * sizeof(double));
   if(temp_A == NULL){
      printf("Memory allocation failed! Aborting function %s.\n", __func__);
      exit(EXIT_FAILURE);
   }

   // Copy in A into temp_A
   for(int rows = 0; rows < m; rows++){
      for(int cols = 0; cols < n; cols++){
         temp_A[rows * (n + m) + cols] = A[rows * n + cols];
      }
   }

   // Copy I into [A I]
   for(int rows = 0; rows < m; rows++){
      for(int cols = 0; cols < m; cols++){
         if(rows == cols){
            temp_A[(n + m) * rows + n + cols] = 1;
         }
         else{
            temp_A[(n + m) * rows + n + cols] = 0;
         }
      }
   }

   // allocate memory for new objective function and modifed minimum.
   double *temp_c = calloc(m + n, sizeof(double));
   if(temp_c == NULL){
      printf("Memory allocation failed! Aborting function %s.\n", __func__);
      exit(EXIT_FAILURE);
   }

   double *y = calloc(n + m, sizeof(double));
   if(y == NULL){
      printf("Memory allocation failed! Aborting function %s. \n", __func__);
      exit(EXIT_FAILURE);
   }

   // copy in values for both temp_c and y.
   for(int i = 0; i < n; i++){
      temp_c[i] = c[i];
      y[i] = x[i];
   }

   // call sf_verify_sol
   bool is_solution_flag = sf_verify_sol(m, m + n, temp_c, temp_A, rhs, y);

   // free allocated memory and clean dangling pointers
   free(temp_A);
   temp_A = NULL;
   free(temp_c);
   temp_c = NULL;
   free(y);
   y = NULL;

   return is_solution_flag;
}
