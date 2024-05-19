#include <cmath>  
#include <fstream>
#include <string>
using namespace std;
#include <stdlib.h>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;
#include <stdint.h>
#include <iostream>
#define profilelog "prof.log"
#define debg "debug.log"
#define outfile "out.txt"
#include <openacc.h>

inline double f(double x)
{
    return std::sin(5*x);
}

double fprime_test(double x)
{
    return 5*std::cos(5*x);
}

// This will save the value of C for viewing, to view only a portion of it, the value of M,N can be changed
void display_1Darr(double* c, int M, int N, string S){
    std::ofstream debugarr(debg, std::ios::app);
    debugarr << S << " " << M << " " << N << endl;
    for (int i=M; i<N;i++){
        debugarr << c[i] << endl;
    }
    debugarr.close();
}

// This will save the value of c, a nx3 array for viewing, to view only a portion of it, the value of M,N can be changed
void display_nx3arr(double* c, int M, int N, string S){
    std::ofstream debugarr(debg, std::ios::app);
    debugarr << S << "from element M " << M << "to element " << N << endl;
    for (int i=M; i<N;i++){
        debugarr << c[3*i] << " " << c[3*i+1] << " " << c[3*i+2] << endl;
    }
    debugarr.close();
}


void saveDoubleArrayToBinary(const double *Array, size_t size, const char *filePath) {
    FILE *file = fopen(filePath, "wb");
    if (file != NULL) {
        fwrite(Array, sizeof(double), size, file);
        fclose(file);
        printf("Float array saved to %s\n", filePath);
    } else {
        fprintf(stderr, "Error: Unable to open file %s for writing.\n", filePath);
    }
}

void pade_solver( int n, double* x)
{
    double xstart=0;
    double xend = 3;
    /*
        n -> number of points to find derivative at

        Solving the differential equation reduces to solving Ax = C, where A is tridiagonal
        THis is done using LU decomposition or Thomas algorithm (TDMA)
    */
    int i;
    double h = (xend - xstart) / (n - 1);
    // creating the 3 diagonals of the matrix a,b,c where b is the main one, in the heap
    //double* a = (double*)acc_malloc((n-1)*sizeof(double));
    //double* b = (double*)acc_malloc(n*sizeof(double));
    //double* c = (double*)acc_malloc(n*sizeof(double));
    //double* v = (double*)acc_malloc(n*sizeof(double));
    
    double* a = (double*)malloc((n-1)*sizeof(double));
    double* b=(double*)malloc(n*sizeof(double));
    double* c=(double*)malloc(n*sizeof(double));
    double* v=(double*)malloc(n*sizeof(double));
    double den;
    #pragma acc data create(a[0:n-1],b[0:n],c[0:n],v[0:n],x[0:n]) copyout(x[0:n]) 
    {
        #pragma acc parallel num_gangs(10) vector_length(32)
        { 
            // Filling the initial and final rows, c has 1 extra 0 element for the algorithm's purpose
            // This is done according to 3rd order accurate pade scheme
            // Filling the remaining rows
            #pragma acc loop
            for (i=1; i<n-1; i++){
                 b[0] = 1;
                b[n-1] = 1;
                a[0] = 2;
                c[n-2] = 2;
                c[n-1] = 0;
                b[i] = 4;
                a[i] = 1;
            }
            #pragma acc loop
            for (i=0; i<n-2; i++){
                c[i] = 1;
            }
            // next store the RHS in vector v
            #pragma acc loop
            for (i = 1; i < n-1; i++)
            {
                v[0] = (-2.5 * f(0) + 2 * f(h) + 0.5 * f(2 * h)) / h;
                v[n - 1] = -(-2.5 * f((n - 1) * h) + 2 * f((n - 2) * h) + 0.5 * f((n - 3) * h)) / h;
                v[i] = 3 * (f(h * (i + 1)) - f(h * (i - 1))) / h;
                c[0] = c[0]/b[0];
                v[0] = c[0]/b[0];
            }
            
            // step1: eliminate the lower diagonal
            
            
            #pragma acc serial
            for (i=1; i<n-1; i++){
                den = b[i]-a[i-1]*c[i-1];
                c[i] = c[i]/den;
                v[i] = (v[i]-a[i-1]*v[i-1])/den;
                x[n-1] = v[n-1];
            }
            // step2: backtrack to fill the values of output array "x", the solution to the linear system
            
            #pragma acc serial
            for (i=n-2; i>-1; i -= 1){
                x[i] = v[i]-c[i]*x[i+1];
            }
        }
    }
    //acc_free(a);
    //acc_free(b);
    //acc_free(c);
    //acc_free(v);
}
int main(int argc, char **argv)
{
    /*
    argv[1] is threadcount
    argv[2] is Number of grid divisions
    */
   
    // Parsing all the command lines arguments
    if (argc > 2)
    {
        cout << "Too many command line arguments format: N" << endl;
        return 1;
    }
    else if (argc < 2)
    {
        cout << "Missing command line arguments format: N" << endl;
        return 1;
    }
    string s = argv[1];
    int N = std::stoi(s);
    // creating array to store output
    double* C=(double*)malloc(N*sizeof(double));

    // beginning timing the program
    auto start_time = Clock::now();
    pade_solver(N, C);
    auto end_time = Clock::now();
    // ending timing the program
    // storing the timing in a log file
    std::ofstream proflog(profilelog, std::ios::app);
    proflog << "N time" << std::endl;
    proflog << N << " " << " " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time    - start_time).count() << endl;
    proflog.close();
    // Writing the output to a txt file for graphing
    std::ofstream outFile(outfile);
    if (outFile.is_open()) {
        for (int i = 0; i < N; ++i) {
            outFile << C[i] << std::endl;
        }
        outFile.close();
        std::cout << "Data written to out.txt successfully." << std::endl;
    } else {
        std::cerr << "Unable to open file out.txt for writing." << std::endl;
    }
    return 0;
}

