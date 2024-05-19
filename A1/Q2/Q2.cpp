#include "Q2.h"

void s_pade_solver(double xstart, double xend, int n, double (*f)(double), vector<double> &x)
{
    /*
        n -> number of points to find derivative at

        Solving the differential equation reduces to solving Ax = C, where A is tridiagonal
        THis is done using LU decomposition or Thomas algorithm (TDMA)
    */
    int i;
    double h = (xend - xstart) / (n - 1);
    // creating the 3 diagonals of the matrix a,b,c where b is the main one, in the heap
    vector<double> a(n-1);
    vector<double> b(n);
    vector<double> c(n);
    // Filling the initial and final rows, c has 1 extra 0 element for the algorithm's purpose
    // This is done according to 3rd order accurate pade scheme
    b[0] = 1;
    b[n-1] = 1;
    a[0] = 2;
    c[n-2] = 2;
    c[n-1] = 0;
    // Filling the remaining rows
    for (i=1; i<n-1; i++){
        b[i] = 4;
        a[i] = 1;
    }
    for (i=0; i<n-2; i++){
        c[i] = 1;
    }
    #ifdef DEBG
    cout << "Tridiagonal a,b,c is made" << endl;
    #endif
    // next store the RHS in vector v
    vector<double> v(n);
    v[0] = (-2.5 * f(0) + 2 * f(h) + 0.5 * f(2 * h)) / h;
    v[n - 1] = -(-2.5 * f((n - 1) * h) + 2 * f((n - 2) * h) + 0.5 * f((n - 3) * h)) / h;
    for (i = 1; i < n-1; i++)
    {
        v[i] = 3 * (f(h * (i + 1)) - f(h * (i - 1))) / h;
    }
    
    // step1: eliminate the lower diagonal
    double den;
    c[0] = c[0]/b[0];
    v[0] = c[0]/b[0];
    for (i=1; i<n-1; i++){
        den = b[i]-a[i-1]*c[i-1];
        c[i] = c[i]/den;
        v[i] = (v[i]-a[i-1]*v[i-1])/den;
    }
    // step2: backtrack to fill the values of output array "x", the solution to the linear system
    x[n-1] = v[n-1];
    for (i=n-2; i>-1; i -= 1){
        x[i] = v[i] -c[i]*x[i+1];
    }
}

void pade_solver( int n, int t, double (*f)(double), vector<double> &x)
{
    // l,d,u form a tridiagonal matrix, Ax = c, l_,d_,_,c_ are its copies
    /*
    double* l = new double(n);
    double* l_ = new double(n);
    double* d = new double(n);
    double* d_ = new double(n);
    double* u = new double(n);
    double* u_ = new double(n);
    double* c = new double(n);
    double* c_ = new double(n);
    */
    vector<double> l(n);
    vector<double> d(n);
    vector<double> u(n);
    vector<double> c(n);
    vector<double> l_(n);
    vector<double> d_(n);
    vector<double> u_(n);
    vector<double> c_(n);

    #pragma omp parallel num_threads(t)
    {
        #pragma omp single
        {
            l[0] = 0; d[0] = 1; u[0] = 2;
            u[n-1] = 0; d[n-1] = 1; l[n-1] = 2;
            c[0] = (n-1)*(-2.5*f(0)+2*f(3.0/(n-1))+0.5*f(6.0/(n-1)));
            c[n-1] = (n-1)*(2.5*f(3)-2*f(3*(1-1.0/(n-1)))-0.5*f(3*(1-2.0/(n-1))));
        }
        #pragma omp for
        for (int i = 1; i<n-1; i++){d[i] = 4;}
        #pragma omp for
        for (int i = 1; i<n-1; i++){l[i] = 1;}
        #pragma omp for
        for (int i = 1; i<n-1; i++){u[i] = 1;}
        #pragma omp for
        for (int i = 1; i<n-1; i++){
            c[i] = (n-1)*(f(3.0*(i+1)/(n-1))-f(3.0*(i-1)/(n-1)));
        }
        for (int k =1; k < n; k *= 2){
            #pragma omp for schedule(dynamic,128) 
            for (int i = 0; i<n; i++){
                if (i>=k && i<n-k){
                    l_[i] = -l[i]*l[i-k]/d[i-k];
                    d_[i] = d[i]-u[i-k]*l[i]/d[i-k]-u[i]*l[i+k]/d[i+k];
                    u_[i] = -u[i]*u[i+k]/d[i+k];
                    c_[i] = c[i]-c[i-k]*l[i]/d[i-k]-c[i+k]*u[i]/d[i+k];
                }
                else if (i>=k && i>=n-k){
                    l_[i] = -l[i]*l[i-k]/d[i-k];
                    d_[i] = d[i]-u[i-k]*l[i]/d[i-k];
                    c_[i] = c[i]-c[i-k]*l[i]/d[i-k];
                }
                else if (i<k && i<n-k){
                    d_[i] = d[i]- u[i]*l[i+k]/d[i+k];
                    u_[i] = -u[i]*u[i+k]/d[i+k];
                    c_[i] = c[i]-c[i+k]*u[i]/d[i+k];
                }
                else{}
            }
            #pragma omp single                  
            {
                l.swap(l_);
                d.swap(d_);
                u.swap(u_);
                c.swap(c_);
            }
        }
        #pragma omp for
        for (int i = 0; i<n; i++){
            x[i] = c[i]/d[i];   
        }
    }
    /*
    delete[] c;
    delete[] c_;
    delete[] l;
    delete[] l_;
    delete[] d;
    delete[] d_;
    delete[] u;
    delete[] u_;
    */
}