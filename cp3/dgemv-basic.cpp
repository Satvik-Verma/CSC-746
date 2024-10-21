const char* dgemv_desc = "Basic implementation of matrix-vector multiply.";

/*
 * This routine performs a dgemv operation
 * Y :=  A * X + Y
 * where A is n-by-n matrix stored in row-major format, and X and Y are n by 1 vectors.
 * On exit, A and X maintain their input values.
 */

void my_dgemv(int n, double* A, double* x, double* y) {
    for (int i = 0; i < n; i++) {
        double temp = 0.0;
        double* A_row = &A[i * n];  // Precompute the start of the row in A
        for (int j = 0; j < n; j++) {
            temp += A_row[j] * x[j];  // Simplified indexing: only use j
        }
        y[i] += temp;
    }
}





