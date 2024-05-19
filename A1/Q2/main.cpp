#include "Q2.h"

int main(int argc, char **argv)
{
    /*
    argv[1] is threadcount
    argv[2] is Number of grid divisions
    */
   
    // Parsing all the command lines arguments
    if (argc > 3)
    {
        cout << "Too many command line arguments format: threads N" << endl;
        return 1;
    }
    else if (argc < 3)
    {
        cout << "Missing command line arguments format: threads N" << endl;
        return 1;
    }
    std::string s(argv[1]);
    int threads = std::stoi(s);
    s = argv[2];
    int N = std::stoi(s);
    // creating array to store output
    vector<double> C(N);

    // beginning timing the program
    double T = omp_get_wtime();
    if (threads==10){
        s_pade_solver(0, 3, N, &f_test, C);
    }
    else{
        pade_solver(N, threads, &f_test, C);
    }
    double T_ = omp_get_wtime();
    // ending timing the program
    // storing the timing in a log file
    std::ofstream proflog(profilelog, std::ios::app);
    proflog << "N threads time" << std::endl;
    proflog << N << " " << threads << " " << T_ - T << endl;
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

double f_test(double x)
{
    return std::sin(5*x);
}

double fprime_test(double x)
{
    return 5*std::cos(5*x);
}
