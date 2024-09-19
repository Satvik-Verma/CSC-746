#include <iostream> // for debugging, can be removed for performance
#include <cstring>  // for memcpy to copy blocks

const char* dgemm_desc = "Blocked dgemm.";

/* This routine performs a blocked/tiling dgemm operation:
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in row-major format.
 * The block size is given by block_size. */
void square_dgemm_blocked(int n, int block_size, double* A, double* B, double* C) 
{
    // Temporary storage for blocks of A, B, and C
    double* A_block = new double[block_size * block_size];
    double* B_block = new double[block_size * block_size];
    double* C_block = new double[block_size * block_size];

    // Iterate over blocks of C
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

    // Free local storage
    delete[] A_block;
    delete[] B_block;
    delete[] C_block;
}