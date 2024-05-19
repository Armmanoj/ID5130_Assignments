#include "Q3.h"
/*
Code design-
The below assumptions are made so as to increase focus on code speed and accuracy rather than on extensive coverage of the topic.
Assumptions and flow of code- Domain is a square, divided into a grid of squares, for solving the PDE. This grid is filled with initial condition, includng boundary,
by the main function, without needing any other script.
A pointer is defined as  a class member that points to this grid.
Then by looping over all non-boundary squares, either using serial (go rw-wise in most cache efficient manner), diagonal or red-black coloring technique (in 3 seperate public functions of the same class Poisson)
, an iteration of the solver is done. Then the max error of the grid is taken, and once this drops below 1 percent, the iterations terminate. The fractional error is 
calculated at each iteration, and compared to see if it is larger than the previous fractional error, as well as  updating the fractional error at beginning of each iteration.
till their ratio x 100 is max percent error is less than 1. 
After executing the code, the result is obtained as the changed values of the grid pointed to. 
This is stored in a txt file (or binary if that is too big,as this is only for graphing, accuracy is less important
Finally, Analysis.py graphs the output, as well as time taken.
Files generated are- an optional debug.log file, an analysis.log file (containing final error, no. of iterations, for given input aguments to main) a profile.log file
(format: delta of grid, threadcount, time taken), 
5 plots (3 3d and 2 time as asked)
A running status of percent error needs to be printed on command line.
Finally, a sehll script is made to automatethis for sir to execte on his machine, and for me to automate it too. Then accuracy and speed up needs to be adjusted 
and understood.


                // checking for the new error in solution
                x = delta*i-1;
                y = delta*j-1;
                err = max(abs(Grid[point]-phi(x,y)), err);
*/
int main(int argc, char** argv)
{
    // reading command line arguments as threads, delta
    if (argc > 4)
    {
        cout << "Too many command line arguments format: threads delta s/rb/d" << endl;
        return 1;
    }
    else if (argc < 4)
    {
        cout << "Missing command line arguments format: threads delta s/rb/d" << endl;
        return 1;
    }
    std::string s(argv[1]);
    int threads = std::stoi(s);
    s = argv[2];
    double delta = std::stod(s);
    s = argv[3];
    // initializing 2 iterators
    int i;
    // seting up the Poisson class and a variable to count iterations
    int iter;
    int N = (int)(2/delta)+1;
    // Variables for timing program
    double t, t_;
    // Creating initial guess of soltion
    double* G = new double[N*N];
    double x,y;
    for (i =0; i<N; i++){
        for (int j=0; j<N; j++){
            G[N*i+j] = 0;
        }
    }

    // Creating RHS function of Poissons equation
    double* q = new double[N*N];
    for (i =0; i<N; i++){
        for (int j =0; j< N; j++){
            x = i*delta-1;
            y = j*delta-1;
            q[i*N+j] = q_(x,y);
        }
    }
    // solving Poisson's equation Serially as well as timing it
    Poisson poisson(delta, N, G, q);
    t = omp_get_wtime();
    if (s=="s"){
        iter = poisson.serial_solve();
    }
    else if (s=="rb"){
        iter = poisson.red_black_solve(threads);
    }
    else if (s=="d"){
        iter = poisson.diagonal_solve(threads);
    }
    else{
        cout << "Give valid solution technique" << endl;
        return 1;
    }
    t_ = omp_get_wtime();
    cout << "delta= " << delta << " Number of iterations = " << iter << endl;
    // storing the timing in a log file
    std::ofstream proflog(profilelog, std::ios::app);
    proflog << "delta threads time" << std::endl;
    proflog << delta << " " << threads << " " << t_-t << endl;
    proflog.close();
    // Saving the solution both as a text file for easy viewing and as a binary file, for preserving the solution accuracy
    saveDoubleArrayToBinary(G, N*N,outbin);
    saveArrayToFile(G, N, outtxt);
    //  Cleaning the pointers and memory allocatiion
    delete [] G;
    delete [] q;
    G=NULL;
    q=NULL;
    poisson.Grid=NULL;
    poisson.q=NULL;
    cout << argv[1] << " " << argv[2] << endl;
    return 0;
}


double q_(double x, double y)
{
    return 2*(2-x*x-y*y);
}
double phi(double x, double y)
{
    return (1-x*x)*(1-y*y);
}