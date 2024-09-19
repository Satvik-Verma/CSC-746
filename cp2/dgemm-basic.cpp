const char* dgemm_desc = "Basic implementation, three-loop dgemm.";

/*
 * This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in row-major format.
 * On exit, A and B maintain their input values.
 */
void square_dgemm(int n, double* A, double* B, double* C) 
{
   // Iterate over rows of C
   for (int i = 0; i < n; ++i) 
   {
      // Iterate over columns of C
      for (int j = 0; j < n; ++j) 
      {
         // Initialize C[i][j] with its current value
         double cij = C[i * n + j]; // C is stored in row-major order
         
         // Perform the dot product of the i-th row of A and j-th column of B
         for (int k = 0; k < n; ++k) 
         {
            cij += A[i * n + k] * B[k * n + j]; // A and B are in row-major order
         }
         
         // Update C[i][j] with the computed value
         C[i * n + j] = cij;
      }
   }
}
