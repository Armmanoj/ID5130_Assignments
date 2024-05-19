/* serial code for Cholesky decomposition */
/* make sure that the init function setups a  */
/* symmetric and positive definite matrix  */
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <cmath>
#define outfile "cholesky.txt"
#define TYPE		float
#define N		10
#define SMALLVALUE	0.001
using namespace std;

void initmult(TYPE mat[][N])
{
  #pragma acc parallel 
  {
  #pragma acc loop
  for (int ii = 0; ii < N; ++ii)
    for (int jj = 0; jj < N && jj < ii; ++jj)
      {	mat[ii][jj] = (ii + jj) / (float)N / N;
	mat[jj][ii] = (ii + jj) / (float)N / N;}

  #pragma acc loop
  for (int ii = 0; ii < N; ++ii)
    mat[ii][ii] = 1.0;
  }
}
			
			
void printMat(TYPE a[][N])
{
  for (int ii = 0; ii < N; ++ii)
    {
      for (int jj = 0; jj < N; ++jj)
	    printf("%.2f ", a[ii][jj]);
      printf("\n");
    }
}

void s_cholesky(TYPE a[][N])
{
  for (int ii = 0; ii < N; ++ii) {
    for (int jj = 0; jj < ii; ++jj) {
      for (int kk = 0; kk < jj; ++kk){
	      a[ii][jj] += -a[ii][kk] * a[jj][kk];
      }
      a[ii][jj] /= (a[jj][jj] > SMALLVALUE ? a[jj][jj] : 1);
      //a[ii][jj] /= a[jj][jj];	// divide by zero.
    }
    for (int kk = 0; kk < ii; ++kk){
      a[ii][ii] += -a[ii][kk] * a[ii][kk];
    }
    a[ii][ii] = sqrt(a[ii][ii]);
  }
}

void cholesky(TYPE a[][N])
{
  #pragma acc data copy(a[0:N][0:N]) 
  // defining a data region like above made the code 10 times faster for N=1000, due to not repeatedly transferring "a"
  {
    #pragma acc parallel num_gangs(N/2) num_workers(1) vector_length(32)
    // putting the parallel region outside the for loop made the code 2.5 times faster, avoid repeated kernel calls
    // the compiler chose a vector length of 128, but 32 is the fastest (2 times better)
    // it is best for a single vector rather than a worker to do the reduction, as only this part has cache efficient memory access
    // if num_gangs is less than N/2, it slows, and same if it is above, the gangs do not play any part in the reduction
    // the time is proportional to num_workers, due to above cache reason
    // the speed up goes from 4.8 times at N=400 to 18 times at N=1000 to 33 times at N=3000
    {
    for (int j=0; j<N;j++){
        TYPE red = 0;
        #pragma acc loop reduction(+:red) 
        for (int k = 0; k<j; k++){
          red += a[j][k]*a[j][k];
        }
        a[j][j] = sqrt(a[j][j]-red);
        #pragma acc loop private(red)
        for (int i=j+1; i<N; i++){
          red = 0;
          #pragma acc loop reduction(+:red) 
          for (int k = 0; k<j; k++){
            red += a[j][k]*a[i][k];
          }
          a[i][j] = a[i][j] - red;
          a[i][j] /= (a[j][j] > SMALLVALUE ? a[j][j] : 1);
        }
      }
    }
  }
}

TYPE errormat(TYPE a[][N], TYPE b[][N]){
  TYPE err = 0;
  for (int i=0; i<N; i++){
    for (int j=0; j<i+1; j++){
      err += fabs(b[i][j]-a[i][j]);
    }
  }
  return err;
}

int main()
{
  TYPE a_s[N][N];
  TYPE a_p[N][N];
  initmult(a_s);
  initmult(a_p);
  //printMat(a_s);
  double T_p;
  double T_p_;
  
  T_p = omp_get_wtime();
  cholesky(a_p);
  T_p_ = omp_get_wtime();
  
  cout << "parallel_time " << T_p_-T_p << "N " << N << endl;
  double T_s = omp_get_wtime();
  s_cholesky(a_s);
  double T_s_ = omp_get_wtime();
  cout << "serial_time " << T_s_-T_s << "N " << N << endl;
  cout << "total_error " << errormat(a_p,a_s) << endl;
  printMat(a_s);
  cout << endl;
  printMat(a_p);
  return 0;
}

