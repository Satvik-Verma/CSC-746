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
   for(int i=0; i<N; i++){
      A[i] = std::rand();
   }
}

int64_t
sum(int64_t N, uint64_t A[])
{
   int64_t sum = 0;
   uint64_t indx = 0;

    while (indx < N) {
        sum += A[indx];
        indx = A[indx];
    }
   printf(" inside sum_indirect perform_sum, N=%lld \n", N);

   return sum;
}

