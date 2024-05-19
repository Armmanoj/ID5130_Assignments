#include "Q3.h"
Poisson::Poisson(double d, int n, double* g, double* Q){
    delta = d;
    N = n;
    Grid = g;
    q = Q;
    for (int i= 0; i<n; i++){
        for (int j=0; j<n; j++){
            q[N*i+j] *= delta*delta;
        }
    }
}


int Poisson::serial_solve() 
{
    // returns the number of iterations run
    int i, j;
    int iterations = 0;
    int point;
    double err = 1;
    while((err > 0.01) && iterations < 1000*N){
        for (i=1; i < N-1; i++){
            for (j=1; j<N-1; j++){
                point = N*i+j;
                Grid[point] = (Grid[point+1]+Grid[point-1]+Grid[point-N]+Grid[point+N]+q[point])*0.25;
            }
        }
        iterations++;
        if (iterations%(5*N) ==  0){
            err=0;
            for (i=1; i < N-1; i++){
                for (j=1; j<N-1; j++){
                    err =  max(err, std::abs(Grid[N*i+j]-(1-(1-(2.0*i)/(N-1))*(1-(2.0*i)/(N-1)))*(1-(1-(2.0*j)/(N-1))*(1-(2.0*j)/(N-1)))));
                }
            } 
        }
    }
    return iterations;
}


int Poisson::red_black_solve(int threads) 
{
    int iterations = 0;
    double err = 1;
    #pragma omp parallel num_threads(threads) 
    {
        int point;
        while((err > 0.01) && iterations < 1000*N){
            
            // red color
            
            #pragma omp for collapse(2) 
            for (int i=1; i < N-1; i+=2){
                for (int j=1; j<N-1; j+=2){
                    point = N*i+j;
                    Grid[point] = (Grid[point+1]+Grid[point-1]+Grid[point-N]+Grid[point+N]+q[point])*0.25;
                }
            }
            #pragma omp for collapse(2) 
            for (int i=2; i < N-1; i+=2){
                for (int j=2; j<N-1; j+=2){
                    point = N*i+j;
                    Grid[point] = (Grid[point+1]+Grid[point-1]+Grid[point-N]+Grid[point+N]+q[point])*0.25;
                }
            }
            // black color
            #pragma omp for collapse(2) 
            for (int i=1; i < N-1; i+=2){
                for (int j=2; j<N-1; j+=2){
                    point = N*i+j;
                    Grid[point] = (Grid[point+1]+Grid[point-1]+Grid[point-N]+Grid[point+N]+q[point])*0.25;
                }
            }
            #pragma omp for collapse(2) schedule(static,64)
            for (int i=2; i < N-1; i+=2){
                for (int j=1; j<N-1; j+=2){
                    point = N*i+j;
                    Grid[point] = (Grid[point+1]+Grid[point-1]+Grid[point-N]+Grid[point+N]+q[point])*0.25;
                }
            }

            #pragma omp single
            {
                iterations++;
                if (iterations%(5*N) ==  0){
                    err=0;
                    for (int i=1; i < N-1; i++){
                        for (int j=1; j<N-1; j++){
                            err =  max(err, std::abs(Grid[N*i+j]-(1-(1-(2.0*i)/(N-1))*(1-(2.0*i)/(N-1)))*(1-(1-(2.0*j)/(N-1))*(1-(2.0*j)/(N-1)))));
                        }
                    } 
                }
            }
        }

    }
    return iterations;
}


int Poisson::diagonal_solve(int threads) 
{
  int iterations = 0;
    double err = 1;
    #pragma omp parallel num_threads(threads) 
    {
        int point;
        while((err > 0.01) && iterations < 500*N){
            // lower and main diagonal
            for (int i=0; i < N-2; i++){
                #pragma omp for 
                for (int j=0; j<=i; j+=1){
                    point = N+1+N*j+i-j;
                    Grid[point] = (Grid[point+1]+Grid[point-1]+Grid[point-N]+Grid[point+N]+q[point])*0.25;
                }
            }
            // upper diagonal
            for (int i=0; i < N-3; i++){
                #pragma omp for
                for (int j=0; j<=i; j+=1){
                    point = N*(N-1)-N*j-i+j-2;
                    Grid[point] = (Grid[point+1]+Grid[point-1]+Grid[point-N]+Grid[point+N]+q[point])*0.25;
                }
            }
            // evaluating error at reglar samples, as well as updating iteration count
            #pragma omp single
            {
                iterations++;
                if (iterations%(5*N) ==  0){
                    err=0;
                    for (int i=1; i < N-1; i++){
                        for (int j=1; j<N-1; j++){
                            err =  max(err, std::abs(Grid[N*i+j]-(1-(1-(2.0*i)/(N-1))*(1-(2.0*i)/(N-1)))*(1-(1-(2.0*j)/(N-1))*(1-(2.0*j)/(N-1)))));
                        }
                    } 
                }
            }
        }
    }
    return iterations;

}