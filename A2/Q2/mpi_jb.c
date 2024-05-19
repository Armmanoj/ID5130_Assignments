#include <mpi.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
# include <math.h>
#define outbin "out.bin"
#define outtxt "out.txt"
#define errbin "err.bin"
#define errtxt "err.txt"

// note that the code needs atleast 2 processors to work

double q_(double x, double y)
{
    return x*x+y*y;
}

void send_boundary(double* G, int M, int N, int rank, int size, MPI_Comm comm){
    // half is either 0 (left half sends data) or 1 (right half sends data)
    if (rank == 0) {
        MPI_Send(G + N * (M - 1) + 1, N - 2, MPI_DOUBLE, rank + 1, 1, comm);
    } 
    else if (rank == size - 1) {
        MPI_Send(G + N+1, N - 2, MPI_DOUBLE, rank - 1, 1, comm);
    }
    else {
        MPI_Send(G + N * M +1, N - 2, MPI_DOUBLE, rank + 1, 1, comm);
        MPI_Send(G + N+1, N - 2, MPI_DOUBLE, rank - 1, 1, comm);
    }
    return;
}

void recv_boundary(double* G, int M, int N, int rank, int size, MPI_Comm comm){
    MPI_Status status;
    if (rank == 0) {
        MPI_Recv(G + N * M + 1, N - 2, MPI_DOUBLE, rank + 1, 1, comm, &status);
    } else if (rank == size - 1) {
        MPI_Recv(G + 1, N - 2, MPI_DOUBLE, rank - 1, 1, comm, &status);     
    } else {
        MPI_Recv(G + N * (M + 1) + 1, N - 2, MPI_DOUBLE, rank + 1, 1, comm, &status);
        MPI_Recv(G + 1, N - 2, MPI_DOUBLE, rank - 1, 1, comm, &status);
    }
    return;
}





void gathererr(double* G, double* Sol, int M, int N, int rank, int size, MPI_Comm Comm){
    // Determine counts and displacements for gather operation
    int* sendcounts = (int*)malloc(size * sizeof(int));
    int* displacements = (int*)malloc(size * sizeof(int));
    int total_count = M * N;

    // Calculate sendcounts based on local value of M on each process
    MPI_Allgather(&M, 1, MPI_INT, sendcounts, 1, MPI_INT, Comm);
    for (int i=0;i<size;i++){
        sendcounts[i]*=N;
    }
    // Calculate displacements
    displacements[0] = 0;
    for (int i = 1; i < size; i++) {
        displacements[i] = displacements[i - 1] + sendcounts[i - 1];
    }

    // Gather data from all processes to processor 0
    MPI_Gatherv(G, total_count, MPI_DOUBLE, Sol, sendcounts, displacements, MPI_DOUBLE, 0, Comm);

    free(sendcounts);
    free(displacements);
    return;
} 

void gatherV(double* G, double* Sol, int M, int N, int rank, int size, MPI_Comm Comm) {
     // Determine counts and displacements for gather operation
    int* sendcounts = (int*)malloc(size * sizeof(int));
    int* displacements = (int*)malloc(size * sizeof(int));
    int total_count = M * N;
    
    // Calculate sendcounts based on local value of M on each process
    MPI_Allgather(&M, 1, MPI_INT, sendcounts, 1, MPI_INT, Comm);
    for (int i=0;i<size;i++){
        sendcounts[i]*=N;
    }

    // Calculate displacements
    displacements[0] = 0;
    for (int i = 1; i < size; i++) {
        displacements[i] = displacements[i - 1] + sendcounts[i - 1];
    }

    // Gather data from all processes to processor 0
    if (rank==0){
        MPI_Gatherv(G , total_count, MPI_DOUBLE, Sol, sendcounts, displacements, MPI_DOUBLE, 0, Comm);
    }
    else{
        MPI_Gatherv(G + N, total_count, MPI_DOUBLE, Sol, sendcounts, displacements, MPI_DOUBLE, 0, Comm);
    }

    free(sendcounts);
    free(displacements);
}

void saveArrayToFile(double* arr, int N, const char* filename) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error opening file for writing.\n");
        return;
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            fprintf(file, "%lf", arr[N * i + j]);
            if (j != N - 1) {
                fprintf(file, " "); // Add space between columns
            }
        }
        fprintf(file, "\n"); // Add newline after each row
    }
    fclose(file);
}


void saveDoubleArrayToBinary(const double *Array, size_t size, const char *filePath) {
    FILE *file = fopen(filePath, "wb");
    if (file != NULL) {
        fwrite(Array, sizeof(double), size, file);
        fclose(file);
    } 
}


          
int Poisson(double* G, double* G1, double* G_err, double* q, MPI_Comm comm, int rank, int size, int M, int N, double* sol) {
    int iter = 0;
    double err = DBL_MAX;
    double errrecv;
    int row_count;
    double* temp;
    if ((rank==0) || (rank==size-1)){
        row_count = M;
    }
    else{
        row_count=M+1;
    }
    MPI_Status status;

    while (err > 0.00000001) {
        for (int i = 1; i < row_count; i++) {
            for (int j = 1; j < N - 1; j++) {
                int idx = i * N + j;
                G1[idx] = (G[(i - 1) * N + j] + G[(i + 1) * N + j] + G[i * N + j - 1] + G[i * N + j + 1] + q[idx-N]) * 0.25;
            }
        }
        for (int i = 1; i < row_count; i++) {
            int idx = i * N + (N - 1);
            G1[idx] = (4 * G1[idx-1] - G1[idx-2]) / 3;
        }
        temp = G1;
        G1 = G;
        G = G1;

        if (rank%2==0){
            send_boundary(G, M, N, rank, size, comm);
        }
        if (rank%2==1){
            recv_boundary(G, M, N, rank, size, comm);
        }
        if (rank%2==1){
            send_boundary(G, M, N, rank, size, comm);
        }
        if (rank%2==0){
            recv_boundary(G, M, N, rank, size, comm);
        }
        if (rank) {
            if ((iter % 100 == 0) && (iter != 0)) {
                for (int i = 0; i < M; i++) {
                    for (int j = 0; j < N; j++) {
                        G_err[i * N + j] = G[(i + 1) * N + j];
                    }
                }
            }
            if ((iter % 100 == 1) && (iter != 1)) {
                err = 0;
                for (int i = 0; i < M; i++) {
                    for (int j = 0; j < N; j++) {
                        err = fmax(err, fabs(G_err[i * N + j] - G[(i + 1) * N + j]));
                    }
                }
            }
        } else {
            if ((iter % 100 == 0) && (iter != 0)) {
                for (int i = 0; i < M; i++) {
                    for (int j = 0; j < N; j++) {
                        G_err[i * N + j] = G[i * N + j];
                    }
                }
            }
            if ((iter % 100 == 1) && (iter != 1)) {
                err = 0;
                for (int i = 0; i < M; i++) {
                    for (int j = 0; j < N; j++) {
                        err = fmax(err, fabs(G_err[i * N + j] - G[i * N + j]));
                    }
                }
            }
        }
        MPI_Allreduce(&err, &errrecv, 1, MPI_DOUBLE, MPI_MAX, comm);
        err = errrecv;
        iter++;
    }
    if (rank) {
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                G_err[i * N + j] = G[(i + 1) * N + j] - G_err[i * N + j];
            }
        }
    } else {
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                G_err[i * N + j] = G[i * N + j] - G_err[i * N + j];
            }
        }
    }
    return iter;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    double I = 0;
    int rank, size;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    // reading command line arguments as delta

    char* s = argv[1];
    char* endptr;
    double delta = strtod(s,&endptr);

    // seting up a variable to count iterations
    int iter;
    // creating array that stores charge distribtion and solution

    int N = (int)(2/delta)+1; // side size of square grid
    int M; // represents the real(without ghost points ) height of array stored in each processor
    if (N%size){
        if (rank<size-1){
            M = N/size+1;
        }
        else{
            M = N-(N/size+1)*(size-1);
        }
    }
    else{
        M = N/size;
    }

// --------------------------------xxx---creating the arrays to be worked with---xxx---------------------------------------------------------
    double* G;
    double* G1;
    if ((rank == 0) || (rank == size-1)){
        G = calloc((M + 1) * N, sizeof(double)); // Flattened array for G with (M + 2) rows
        G1 = calloc((M + 1) * N, sizeof(double));
    }
    else{
        G = calloc((M + 2) * N, sizeof(double)); // Flattened array for G with (M + 2) rows
        G1 = calloc((M + 2) * N, sizeof(double));
    }
    double* G_err = calloc(M * N, sizeof(double)); // Flattened array for G_err with M rows

    if (rank == 0) {
        for (int i = 1; i < M + 1; i++) {
            G[i * N] = sin(2 * M_PI * (i * delta-1));
            G1[i * N] = sin(2 * M_PI * (i * delta-1));
        }
    } else if (rank == size - 1) {
        if (N % size) {
            for (int i = 0; i < M; i++) {
                G[i * N] = sin(2 * M_PI * (i + (N / size + 1) * rank - 1) * delta);
                G1[i * N] = sin(2 * M_PI * (i + (N / size + 1) * rank - 1) * delta);
            }
        } else {
            for (int i = 1; i < M + 1; i++) {
                G[i * N] = sin(2 * M_PI * (i + (N / size) * rank - 1) * delta);
                G1[i * N] = sin(2 * M_PI * (i + (N / size) * rank - 1) * delta);
            }
        }
    } else {
        for (int i = 1; i < M + 1; i++) {
            G[i * N] = sin(2 * M_PI * (i + M * rank - 1) * delta);
            G1[i * N] = sin(2 * M_PI * (i + M * rank - 1) * delta);
        }
    }

    double* q = calloc(M * N, sizeof(double)); // Flattened array for charge distribution, no need for ghost rows here
    if (rank == 0) {
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                q[i * N + j] = q_(j * delta-1, i * delta-1)*delta*delta;
            }
        }
    } else if (rank == size - 1) {
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                if (N % size) {
                    q[i * N + j] = q_(j * delta-1, (i  + rank * (N / size + 1))*delta-1)*delta*delta;
                } else {
                    q[i * N + j] = q_(j * delta-1, (i + N - N / size)*delta-1)*delta*delta;
                }
            }
        }
    } else {
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                q[i * N + j] = q_(j * delta-1, (i + M * rank)*delta-1)*delta*delta;
            }
        }
    }

    // creating a buffer array to hold the full solution grid for writing to output file
    double* Sol = NULL;
    double* Err = NULL;
    if (rank==0){
        Sol = (double*)malloc(N*N*sizeof(double));
        Err = (double*)malloc(N*N*sizeof(double));
    }
//---------------------------------------------xxx--------------------------------xxx---------------------------------------------------------


    // solving Poisson's equation as well as timing it
    double start = MPI_Wtime();
    iter = Poisson(G,G1,G_err,q,MPI_COMM_WORLD,rank,size,M,N,Sol);
    double end = MPI_Wtime();\
    // saving log data
    double time = end-start;
    double time_recv;
    MPI_Reduce( &time , &time_recv ,1 ,MPI_DOUBLE , MPI_MAX , 0 , MPI_COMM_WORLD);
    if (rank == 0){
        printf("%f\n",time_recv);//delta= %f Number of iterations = %d processors = %d time = %f",delta,iter,size,time_recv) ;
    }

    // Gathering the local solution arrays together
    gathererr(G_err,Err, M, N,rank, size, MPI_COMM_WORLD);
    MPI_Barrier( MPI_COMM_WORLD);
    gatherV(G, Sol, M, N,rank, size, MPI_COMM_WORLD);
    
    if (rank==0){
        // Saving the solution both as a text file for easy viewing and as a binary file, for preserving the solution accuracy
        saveDoubleArrayToBinary(Sol, N*N,argv[2]);
        saveArrayToFile(Sol, N, outtxt);
        saveDoubleArrayToBinary(Err, N*N,errbin);
        saveArrayToFile(Err, N, errtxt);
    }
    //  Cleaning the pointers and memory allocatiion
    free(G);
    free(q);
    free(G_err);
    if (rank==0){
        free(Sol);
        Sol = NULL;
    }
    G=NULL;
    G_err=NULL;
    q=NULL;
    MPI_Finalize();
    return 0;
}


