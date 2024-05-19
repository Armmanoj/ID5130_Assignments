#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
using namespace std;

// defining files used
#define profilelog "profile.log"
#define debg "debug.log"
#define outfile "out.txt"

// Primary functions
void s_pade_solver(double xstart, double xend, int n, double (*f)(double), vector<double> &x);
void pade_solver(int n, int t, double (*f)(double), vector<double> &C);

// Testing functions
double f_test(double x);
double fprime_test(double x);

// Utilities
typedef struct{
    int start, stop;
}thread_range;
thread_range t_range(int N, int threads, int id);

void display_1Darr(double* c, int M, int N, string S); // saving 1D arrays to txt file 
void display_nx3arr(double* c, int M, int N, string S); // saving nx3 matrices to txt file
void saveDoubleArrayToBinary(const double *Array, size_t size, const char *filePath); //saving double arrays to binary file

