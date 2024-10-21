#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include "likwid-stuff.h"

const char* dgemm_desc = "Blocked dgemm, OpenMP-enabled";


/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in row-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm_blocked(int n, int block_size, double* A, double* B, double* C) 
{
   // insert your code here: implementation of blocked matrix multiply with copy optimization and OpenMP parallelism enabled

   // be sure to include LIKWID_MARKER_START(MY_MARKER_REGION_NAME) inside the block of parallel code,
   // but before your matrix multiply code, and then include LIKWID_MARKER_STOP(MY_MARKER_REGION_NAME)
   // after the matrix multiply code but before the end of the parallel code block.

   std::cout << "Insert your blocked matrix multiply with copy optimization, openmp-parallel edition here " << std::endl;
    // Temporary storage for blocks of A, B, and C
    double* A_block = new double[block_size * block_size];
    double* B_block = new double[block_size * block_size];
    double* C_block = new double[block_size * block_size];

    // Start LIKWID marker region
    LIKWID_MARKER_START("blocked");

    // Parallelize the outer loops using OpenMP
    #pragma omp parallel for collapse(2) private(A_block, B_block, C_block)
    for (int i = 0; i < n; i += block_size)
    {
        for (int j = 0; j < n; j += block_size)
        {
            // Determine the actual block size for the current block
            int M = std::min(block_size, n - i);
            int N = std::min(block_size, n - j);

            // Initialize C_block with the current values of C
            for (int ii = 0; ii < M; ++ii)
                for (int jj = 0; jj < N; ++jj)
                    C_block[ii * block_size + jj] = C[(i + ii) * n + (j + jj)];

            // Iterate over blocks of A and B
            for (int k = 0; k < n; k += block_size)
            {
                // Determine the actual block size for the current block
                int K = std::min(block_size, n - k);

                // Load block of A into A_block (row-major)
                for (int ii = 0; ii < M; ++ii)
                    for (int kk = 0; kk < K; ++kk)
                        A_block[ii * block_size + kk] = A[(i + ii) * n + (k + kk)];

                // Load block of B into B_block (row-major)
                for (int kk = 0; kk < K; ++kk)
                    for (int jj = 0; jj < N; ++jj)
                        B_block[kk * block_size + jj] = B[(k + kk) * n + (j + jj)];

                // Perform block matrix multiplication: C_block += A_block * B_block
                for (int ii = 0; ii < M; ++ii) 
                {
                    for (int jj = 0; jj < N; ++jj) 
                    {
                        for (int kk = 0; kk < K; ++kk) 
                        {
                            // row-major indexing
                            C_block[ii * block_size + jj] += A_block[ii * block_size + kk] * B_block[kk * block_size + jj];
                        }
                    }
                }
            }

            // Copy the updated C_block back into the C matrix (row-major)
            for (int ii = 0; ii < M; ++ii)
                for (int jj = 0; jj < N; ++jj)
                    C[(i + ii) * n + (j + jj)] = C_block[ii * block_size + jj];
        }
    }

    // Stop LIKWID marker region
    LIKWID_MARKER_STOP("blocked");

    // Free local storage
    delete[] A_block;
    delete[] B_block;
    delete[] C_block;
}
