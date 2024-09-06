#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
#include <string.h>



void 
setup(int64_t N, uint64_t A[])
{
   printf(" inside sum_indirect problem_setup, N=%lld \n", N);
   std::random_device rd;
   std::mt19937 gen(rd());
   std::uniform_int_distribution<uint64_t> dis(0, N-1);

    // Fill array A with random values
   for (int64_t i = 0; i < N; ++i) {
      A[i] = dis(gen);
   }
}

int64_t
sum(int64_t N, uint64_t A[])
{
    int64_t sum = 0;

    // Compute the sum accumulating A[indx] where indx = A[indx]
    for (int64_t i = 0; i < N; ++i) {
        sum += A[A[i]];
    }

   printf(" inside sum_indirect perform_sum, N=%lld \n", N);

   return sum;
}

