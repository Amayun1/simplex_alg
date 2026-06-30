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
double sf_simplex(lapack_int m, lapack_int n, const double *c, const double *A, const double *rhs, volatile double *bfs);
void sf_simplex_phase_one(lapack_int m, lapack_int n, const double *A, const double *rhs, volatile double *bfs);
double ez_simplex(lapack_int m, lapack_int n, const double *c, const double *A, const double *rhs, volatile double *x);
void print_float_matrix(int m, int n, double *matrix);
void print_bool_matrix(int m, int n, bool *matrix);
void print_int_matrix(int m, int n, int *matrix);
bool is_basic(int m, int n, const double *matrix, const double *rhs, const double *bfs);
