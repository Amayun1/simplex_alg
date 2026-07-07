/*
Header file for simplex algorithm library written in C.
*/

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <lapacke.h>
#include <cblas.h>

// Header for tolower command for user input in the warning.
#include <ctype.h>

// Header for INFINITY macro.
#include <math.h>

// define the maximum number of iterations the algorithm is allowed to do.
#ifndef MAX_SIMPLEX_ITERS
#define MAX_SIMPLEX_ITERS 10000
#endif

// Function Prototypes here
double sf_simplex(lapack_int m, lapack_int n, double *c, double *A, double *rhs, double *bfs, bool *BFS_indices);
void sf_simplex_phase_one(lapack_int m, lapack_int n, double *A, double *rhs, double *bfs, bool *BFS_indices);
double ez_simplex(lapack_int m, lapack_int n, double *c, double *A, double *rhs, double *x);
bool is_basic(int m, int n, double *A, double *rhs, double *bfs);
bool sf_verify_sol(int m, int n, double *c, double *A, double *rhs, double *x);
