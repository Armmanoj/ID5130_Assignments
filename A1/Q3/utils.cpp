#include "Q3.h"

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


void saveDoubleArrayToBinary(const double *Array, size_t size, const char *filePath) {
    FILE *file = fopen(filePath, "wb");
    if (file != NULL) {
        fwrite(Array, sizeof(double), size, file);
        fclose(file);
        cout << "Float array saved to " << filePath << endl;
    }
    else {
        cout << "Error: Unable to open file " << filePath << " for writing" << endl;
    }
}

void saveArrayToFile(double* arr, int N, const char* filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file for writing." << std::endl;
        return;
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            file << arr[N*i+j];
            if (j != N - 1) {
                file << " "; // Add space between columns
            }
        }
        file << std::endl; // Add newline after each row
    }
    cout << "Float array saved to " << filename << endl;
    file.close();
}