/*

Functions for printing some matrices.

*/

void print_float_matrix(int m, int n, double *matrix){
   //prints matrix in row major order with m rows and n columns.
   for(int row = 0; row < m; row++){
      for(int column = 0; column < n; column++){
         printf("%10.6f", matrix[row * n + column]);
      }
      printf("\n");
   }
   printf("\n");
   return;
}

void print_bool_matrix(int m, int n, bool *matrix){
   //prints matrix in row major order with m rows and n columns.
   for(int row = 0; row < m; row++){
      for(int column = 0; column < n; column++){
         printf("%10d", matrix[row * n + column]);
      }
      printf("\n");
   }
   printf("\n");
   return;
}

void print_int_matrix(int m, int n, int *matrix){
   //prints matrix in row major order with m rows and n columns.
   for(int row = 0; row < m; row++){
      for(int column = 0; column < n; column++){
         printf("%10d", matrix[row * n + column]);
      }
      printf("\n");
   }
   printf("\n");
   return;
}

void print_long_matrix(int m, int n, long *matrix){
   //prints matrix in row major order with m rows and n columns.
   for(int row = 0; row < m; row++){
      for(int column = 0; column < n; column++){
         printf("%10ld", matrix[row * n + column]);
      }
      printf("\n");
   }
   printf("\n");
   return;
}

void print_longlong_matrix(int m, int n, long long *matrix){
   //prints matrix in row major order with m rows and n columns.
   for(int row = 0; row < m; row++){
      for(int column = 0; column < n; column++){
         printf("%10lld", matrix[row * n + column]);
      }
      printf("\n");
   }
   printf("\n");
   return;
}

// combine the above into one boiler-plate print matrix function. Only available in C11 or higher due to _Generic
#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L

#define print_matrix(m, n, matrix) _Generic((matrix),	\
   int*		: print_int_matrix,		\
   long*	: print_long_matrix,		\
   long long*	: print_longlong_matrix,	\
   float*	: print_float_matrix,		\
   double*	: print_float_matrix,		\
   bool*	: print_bool_matrix,		\
   default	: printf("Error! Tried to print a matrix of unknown type. See print_matrix.c for compatable types.\n") \
)( (m), (n), (matrix) )

#endif

