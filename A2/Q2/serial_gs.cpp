#include <vector>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <chrono>
using namespace std;
#define outbin "out.bin"
#define errbin "err.bin"
#define outtxt "out.txt"
#define errtxt "err.txt"

inline double q_(double x, double y)
{
    return x*x+y*y;
}

void saveArrayToFile(const vector<
vector<double>>& arr, const char* filename) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error opening file for writing.\n");
        return;
    }

    for (const auto& row : arr) {
        for (size_t j = 0; j < row.size(); ++j) {
            fprintf(file, "%lf", row[j]);
            if (j != row.size() - 1) {
                fprintf(file, " "); // Add space between columns
            }
        }
        fprintf(file, "\n"); // Add newline after each row
    }
    fclose(file);
}

void saveDoubleArrayToBinary(const std::vector<std::vector<double>>& Array, const char* filePath) {
    FILE* file = fopen(filePath, "wb");
    if (file != NULL) {
        for (const auto& row : Array) {
            fwrite(row.data(), sizeof(double), row.size(), file);
        }
        fclose(file);
    } else {
        std::cout << "Error: Unable to open file " << filePath << " for writing" << std::endl;
    }
}

int Poisson(vector<vector<double>>& G, vector<vector<double>>& q,int N, double delta,vector<vector<double>>& G_err){
    int iter = 0;
    double err = MAXFLOAT;
    while (err>0.00000001){
        for (int i=1; i<N-1; i++){
            for (int j=1; j<N-1; j++){
                G[i][j] = (G[i-1][j]+G[i+1][j]+G[i][j-1]+G[i][j+1]+q[i][j])*0.25;
            }
            G[i][N-1] = (4*G[i][N-2]-G[i][N-3])/3;
        }
        iter++;
        if ((iter%100 == 0) && (iter!=0)){
            for (int i=0; i<N; i++){
                for (int j=0; j<N; j++){
                    G_err[i][j] = G[i][j];
                }
            }
        }
        if ((iter%100 == 1) && (iter!=1)){
            err = 0;
            for (int i=0; i<N; i++){
                for (int j=0; j<N; j++){
                    err = max(err,abs(G_err[i][j] - G[i][j]));
                }
            }
        }
    }
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            G_err[i][j] = G[i][j] - G_err[i][j];
        }
    }
    return iter;
}

int main(int argc, char** argv){
    /*
     main first initializes G and q,= with a given delta, N. Then runs Poissons, runs timing and stores it     */
    string s = argv[1];
    double delta = std::stod(s);
    int N = 2/delta+1;
    vector<vector<double>> q(N, vector<double>(N, 0));
    for (int i =0; i<N; i++){
        double x,y;
        for (int j =0; j< N; j++){
            x = i*delta-1;
            y = j*delta-1;
            q[i][j] = delta*delta*q_(x,y);
        }
    }
    vector<vector<double>> G(N, vector<double>(N, 0));
    vector<vector<double>> G_err(N, vector<double>(N, 0));
    for (int j=0; j<N;j++){
        G[j][0] = sin(2*M_PI*(j*delta-1));
    }
    auto start = chrono::high_resolution_clock::now();
    int iter = Poisson(G,q,N,delta,G_err);
    auto end = std::chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;
    cout << "delta= " << delta << " Number of iterations = " << iter <<" processors = " << 1 << " time = " << duration.count()  << std::endl;
    saveDoubleArrayToBinary(G, outbin);
    saveDoubleArrayToBinary(G_err, errbin);
    saveArrayToFile(G, outtxt);
    saveArrayToFile(G, errtxt);
    return 0;
}
 