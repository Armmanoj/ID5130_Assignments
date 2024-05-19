#include <stdio.h>
#include <omp.h>

int main() {
  #ifdef _OPENMP
    printf("OpenMP version: %d\n", _OPENMP);
  #else
    printf("No OpenMP support\n");
  #endif
  return 0;
}