const char* dgemm_desc = "Basic implementation, three-loop dgemm.";

/*
 * This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in row-major format.
 * On exit, A and B maintain their input values.
 */
void square_dgemm(int n, double* A, double* B, double* C) 
{
   // insert your code here: implementation of basic matrix multiple where A and B are in row-major format
   for (int i = 0; i < n; ++i)
   {
      for (int j = 0; j < n; ++j)
      {
         for (int k = 0; k < n; ++k)
         {
            C[i+j*n] += A[i+k*n] * B[k+j*n];
         }
      }
   }
}
