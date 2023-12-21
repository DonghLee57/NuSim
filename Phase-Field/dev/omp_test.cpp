#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <ctime>
#include <omp.h>
using namespace std;

int main()
{
    double start, end;
    double M[100][100] = {1};
    int SUM = 0;
    // Evolve
    start = omp_get_wtime();
    //int N_THR = omp_get_num_threads();
    //cout << to_string(omp_get_num_threads()) << " threads\n";
    //cout << to_string(omp_get_thread_num()) << " thread id\n";
    #pragma omp parallel for
    for (int k=0; k <= 100000; k++) {
       //#pragma omp parallel for
       for (int x=0; x<100; x++) {
            for (int y=0; y<100; y++) {
                SUM = SUM + M[x][y];
            }
       }
       //#pragma omp critical
       //if (k == 0) cout << to_string(omp_get_num_threads()) << " threads\n";
    } // MAXIT loop
    end = omp_get_wtime();
    cout << SUM << endl;
    cout << (double)(end - start) << endl;
    return 0; 
}
