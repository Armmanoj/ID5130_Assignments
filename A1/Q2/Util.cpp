#include "Q2.h"

thread_range t_range(
    int N, int threads, int id){
    int start, stop;
    if (N%threads == 0){
        start = id*(N/threads);
        stop = std::min(N, (id+1)*(N/threads)); 
    }
    else{
        start = id*(N/(threads-1));          
        stop = std::min(N, (id+1)*(N/(threads-1)));
    }
    return thread_range{start,stop};
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