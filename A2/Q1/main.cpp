#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#define dx 0.00002
#define dt 0.00002
#define c 1
#define L 2
#define T 1
#define fileout "out.txt"
using namespace std;
// boundaries are 0 to L in x, and 0 to T in t

vector<double> upwind(){
    int N = L/dx+1; 
    int M = T/dt+1;
    vector<double> u(N*M,0);
    // now setting initial conditions
    for (int i=0; i<0.5/dx;i++){
        u[i] = sin(4*M_PI*i*dx);
    }
    double sol_speed = c*dt/dx;
    for(int t=0;t<M-1;t+=1){
        for(int i=1; i<N-1;i++){
            u[i+N*(t+1)] = u[i+N*t]-sol_speed*(u[i+N*t]-u[i+N*t-1]);
        }
    }
    return u;
}

vector<double> QUICK(){
    int N = L/dx+1;
    int M = T/dt+1;
    vector<double> u(N*M,0);
    // now setting initial conditions
    for (int i=0; i<0.5/dx;i++){
        u[i] = sin(4*M_PI*i*dx);
    }
    double sol_speed = c*dt/dx;
    for(int t=0;t<M-1;t+=1){
        u[1+N*(t+1)] = u[1+N*t]-sol_speed*(u[1+N*t]-u[N*t]);
        for(int i=2; i<N-1;i++){
            u[i+N*(t+1)] = u[i+N*t]-sol_speed*((u[i+N*t]+u[i+N*t+1])*0.375-0.875*u[i+N*t-1]+0.125*u[i+N*t-2]);
        }
    }
    return u;
}


void store_data(const std::string& filename, const std::vector<double>& u, int N, int M) {
    // Open the file for writing
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
    }

    // Write the data to the file
    for (int j = 0; j<M; j++) {
        if ((j==0) || (j==M/2) || (j==M-1)){
            for (int i = 0; i < N; ++i) {
            file << u[j*N+i]; // Write element
            if (i<N-1) {
                file << " "; // Separate elements by space, except for the last one in each row
            }
        }
        file << "\n"; // Newline after each row
        } 
    }

    // Close the file
    file.close();
}

int main(int argc, char** argv){
    int N = L/dx+1;
    int M = T/dt+1;   
    if (argv[1][0] == 'u'){
        store_data(fileout,upwind(),N,M);
    }
    if (argv[1][0] == 'q'){
        store_data(fileout,QUICK(),N,M);
    }
    return 0;
}