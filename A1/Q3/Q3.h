#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <fstream>
#include <string>
#include <iostream>
using namespace std;

// files defined
#define outbin "out.bin"
#define outtxt "out.txt"
#define profilelog "prof.log"

// Test functions
double q_(double x, double y);
double phi(double x, double y);

// Main solver class POisson
class Poisson{
    public:
        int N; // side of square, centered at origin
        double* Grid; // pointer to grid
        double* q; // pointer to grid that stores value of the RHS of Poissons equation at each point
        double delta;
        int point;
        // constructor
        Poisson(double d, int n, double* g, double* Q);
        //solvers
        int serial_solve();
        int red_black_solve(int threads);
        int diagonal_solve(int threads);
};

// Utilities
typedef struct{
    int start, stop;
}thread_range;
thread_range t_range(int N, int threads, int id);

void saveDoubleArrayToBinary(const double *Array, size_t size, const char *filePath) ;

void saveArrayToFile(double* arr, int N, const char* filename);