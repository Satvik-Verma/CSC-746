const char* dgemv_desc = "Basic implementation of matrix-vector multiply.";

/*
 * This routine performs a dgemv operation
 * Y :=  A * X + Y
 * where A is n-by-n matrix stored in row-major format, and X and Y are n by 1 vectors.
 * On exit, A and X maintain their input values.
 */
void my_dgemv(int N, double *A, double *x, double *y) {
    for (int i = 0; i < N; ++i) {
        double temp = 0.0; // Accumulate the dot product
        double *A_row = &A[i * N]; // Pointer to the start of the current row of A
        for (int j = 0; j < N; ++j) {
            temp += A_row[j] * x[j]; // Access elements of A using a simple stride
        }
        y[i] += temp;
    }
}



