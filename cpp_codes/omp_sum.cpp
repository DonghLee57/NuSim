#include <iostream>
#include <omp.h>

int main(){
    int Nx = 100, Ny = 100;
    // Simple summation
    int count = 0;
    for (int x = 1; x < Nx + 1; ++x) 
    for (int y = 1; y < Ny + 1; ++y) 
        count++;
    printf("Simple computation\n");
    printf("Count is %d\n", count);

    // OpenMP summation
    count = 0;
    #pragma omp parallel for
    for (int x = 1; x < Nx + 1; ++x) 
    for (int y = 1; y < Ny + 1; ++y) 
        count++;
    printf("OpenMP without reduction\n");
    printf("Count is %d\n", count);

    count = 0;
    #pragma omp parallel for reduction(+ : count)
    for (int x = 1; x < Nx + 1; ++x) 
    for (int y = 1; y < Ny + 1; ++y) 
        count++;
    printf("OpenMP with reduction\n");
    printf("Count is %d\n", count);

    return 0;
}

/* - Example Output -
Simple computation
Count is 100
OpenMP without reduction
Count is 70
OpenMP with reduction
Count is 100
*/
