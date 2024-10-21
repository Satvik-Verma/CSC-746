const char* dgemv_desc = "Vectorized implementation of matrix-vector multiply.";

/*
 * This routine performs a dgemv operation
 * Y := A * X + Y
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



// /global/homes/s/satvik/CSC-746/cp3/dgemv-vectorized.cpp:10:23: missed: couldn't vectorize loop
// /global/homes/s/satvik/CSC-746/cp3/dgemv-vectorized.cpp:10:23: missed: not vectorized: loop nest containing two or more consecutive inner loops cann>
// /global/homes/s/satvik/CSC-746/cp3/dgemv-vectorized.cpp:21:18: optimized: loop vectorized using 16 byte vectors
// /global/homes/s/satvik/CSC-746/cp3/dgemv-vectorized.cpp:21:18: optimized:  loop versioned for vectorization because of possible aliasing
// /global/homes/s/satvik/CSC-746/cp3/dgemv-vectorized.cpp:14:18: optimized: loop vectorized using 32 byte vectors
// /global/homes/s/satvik/CSC-746/cp3/dgemv-vectorized.cpp:14:18: optimized:  loop versioned for vectorization because of possible aliasing
// /global/homes/s/satvik/CSC-746/cp3/dgemv-vectorized.cpp:14:18: optimized: loop vectorized using 16 byte vectors
// /global/homes/s/satvik/CSC-746/cp3/dgemv-vectorized.cpp:9:6: note: vectorized 2 loops in function.
// /global/homes/s/satvik/CSC-746/cp3/dgemv-vectorized.cpp:9:6: note: ***** Analysis failed with vector mode V4DF
// /global/homes/s/satvik/CSC-746/cp3/dgemv-vectorized.cpp:9:6: note: ***** The result for vector mode V32QI would be the same
// /global/homes/s/satvik/CSC-746/cp3/dgemv-vectorized.cpp:9:6: note: ***** Re-trying analysis with vector mode V16QI
// /global/homes/s/satvik/CSC-746/cp3/dgemv-vectorized.cpp:9:6: note: ***** Analysis failed with vector mode V16QI
// /global/homes/s/satvik/CSC-746/cp3/dgemv-vectorized.cpp:9:6: note: ***** Re-trying analysis with vector mode V8QI
// /global/homes/s/satvik/CSC-746/cp3/dgemv-vectorized.cpp:9:6: note: ***** Analysis failed with vector mode V8QI
// /global/homes/s/satvik/CSC-746/cp3/dgemv-vectorized.cpp:9:6: note: ***** Re-trying analysis with vector mode V4QI
// /global/homes/s/satvik/CSC-746/cp3/dgemv-vectorized.cpp:9:6: note: ***** Analysis failed with vector mode V4QI
// void my_dgemv(int N, double *A, double *x, double *y) {
//     for (int j = 0; j < N; ++j) {
//         double x_j = x[j];
//         int i = 0;

//         // Unroll the loop to process 4 elements at a time
//         for (; i <= N - 4; i += 4) {
//             double temp0 = y[i];
//             double temp1 = y[i + 1];
//             double temp2 = y[i + 2];
//             double temp3 = y[i + 3];

//             temp0 += A[i * N + j] * x_j;
//             temp1 += A[(i + 1) * N + j] * x_j;
//             temp2 += A[(i + 2) * N + j] * x_j;
//             temp3 += A[(i + 3) * N + j] * x_j;

//             y[i]     = temp0;
//             y[i + 1] = temp1;
//             y[i + 2] = temp2;
//             y[i + 3] = temp3;
//         }

//         // Process any remaining elements
//         for (; i < N; ++i) {
//             y[i] += A[i * N + j] * x_j;
//         }
//     }
// }