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
    double A_block[block_size * block_size];
    double B_block[block_size * block_size];
    double C_block[block_size * block_size];

    // Iterate over blocks of C
    for (int i = 0; i < n; i += block_size)
    {
        for (int j = 0; j < n; j += block_size)
        {
            // Copy the current block of C into C_block
            for (int ii = 0; ii < block_size; ++ii)
                memcpy(&C_block[ii * block_size], &C[(i + ii) * n + j], block_size * sizeof(double));

            // Iterate over blocks of A and B
            for (int k = 0; k < n; k += block_size)
            {
                // Load block of A into A_block (row-major)
                for (int ii = 0; ii < block_size; ++ii)
                    memcpy(&A_block[ii * block_size], &A[(i + ii) * n + k], block_size * sizeof(double));

                // Load block of B into B_block (row-major)
                for (int kk = 0; kk < block_size; ++kk)
                    memcpy(&B_block[kk * block_size], &B[(k + kk) * n + j], block_size * sizeof(double));

                // Perform block matrix multiplication: C_block += A_block * B_block
                for (int ii = 0; ii < block_size; ++ii) 
                {
                    for (int jj = 0; jj < block_size; ++jj) 
                    {
                        double sum = C_block[ii * block_size + jj]; // Load C_block element
                        for (int kk = 0; kk < block_size; ++kk) 
                        {
                            // row-major indexing
                            sum += A_block[ii * block_size + kk] * B_block[kk * block_size + jj];
                        }
                        C_block[ii * block_size + jj] = sum; // Store the result
                    }
                }
            }

            // Copy the updated C_block back into the C matrix (row-major)
            for (int ii = 0; ii < block_size; ++ii)
                memcpy(&C[(i + ii) * n + j], &C_block[ii * block_size], block_size * sizeof(double));
        }
    }
}