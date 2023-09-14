#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <random>
using namespace std;

// Discretization parameters
const int MAXIT = 10;   // Iteration number
const int NSTEP = 1;    // Step size for save
const int Ngr=    1;    // The number of grains
const int Nx = 10;      // Discretization x-direction
const int Ny = 10;      // Discretization y-direction
double dx = 0.1;        // Grid spacing in x-direction [Ang]
double dy = 0.1;        // Grid spacing in y-direction [Ang]
double dt = 0.001;      // Size of time step [s]
double Mob=  10.0;             // Mobility
double grad_coeff = 0.5;       // gradient coefficient
double A = 1, B = 1, Gam = 1;
const int BCells = 1; 		     // Number of boundary cells

// Phase-field storage
const int Mx = Nx + 2 * BCells;
const int My = Ny + 2 * BCells;

double ETAS     [Ngr][Mx][My] = {0};
double ETAS_new [Ngr][Mx][My] = {0};
double df       [Mx][My] = {0};
double Lap_ETAS [Mx][My] = {0};

void SAVE(int step)
{
    ofstream myfile("PhaseField_" + to_string(step) + ".dat");
    for (int x=1; x<Nx+1; x++) {
        for (int y=1; y<Ny+1; y++) {
        myfile << double(ETAS[Ngr][x][y]) << "  ";
        }
        myfile << "\n";
    }
    myfile.close();
}

int main()
{
    // Initiallization
    for (int n=1; n <= Ngr; n++) {
        int xc = rand() % (Nx+1) + 1;
        int yc = rand() % (Ny+1) + 1;
        for (int x=1; x<Nx+1; x++) {
            for (int y=1; y<Ny+1; y++) {
                if (sqrt(pow(xc-x,2) + pow(yc-y,2)) < 2)
                {
                    ETAS[n][x][y] = 1.0;
                }
            }
        } 
    }

    // Evolve
    for (int k=1; k <= MAXIT; k++) {
//        for (int n=1; n <= Ngr; n++) {
//        }
        SAVE(k);
    }
    return 0; 
}
