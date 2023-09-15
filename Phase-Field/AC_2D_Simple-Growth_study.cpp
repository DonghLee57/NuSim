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

// Functions
double lap_2D(double Grid, int x, int y, double dx, double dy, int Nx, int Ny);
void SAVE(int step, int Ngr);

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
        // Save data
        if (k % NSTEP == 0){
            for (int n=0; n < Ngr; n++)
            SAVE(k, n);
        }
        // Calculate
        for (int N=0; N < Ngr; N++){
            Lap_ETAS[Mx][My] = {0};
            df[Mx][My] = {0};
            
            for (int x=1; x<Nx+1; x++) {
            for (int y=1; y<Ny+1; y++) {
                Lap_ETAS[x][y] = lap_2D(ETAS[N], x, y, dx, dy, Nx, Ny);
                int SUM_eta_sq = 0;
                for (int ngr=0; ngr < Ngr; ngr++){
                    if (N != ngr){
                    SUM_eta_sq = SUM_eta_sq + ETAS[ngr][x][y]**2);
                    }
                    df[x][y] = -A*ETAS[N][x][y] + B*ETAS[N][x][y]**3 + 2*ETAS[N][x][y]*SUM_eta_sq;
                }
            }
            }
            for (int x=1; x<Nx+1; x++) {
            for (int y=1; y<Ny+1; y++) {
                ETAS_new[N][x][y] = ETAS[N][x][y] - dt*Mob*(df[x][y] - grad_coeff*Lap_ETAS[x][y]);
                }
            }
        // check result
        ETAS = ETAS_new;
        }
    }
    return 0; 
}

double lap_2D(double Grid, int x, int y, double dx, double dy, int Nx, int Ny)
{
    int idx0 = x - 1;
    int idxN = x + 1;
    int idy0 = y - 1;
    int idyN = y + 1;
    res = (Grid[idx0][y] + Grid[idxN][y] - 2*Grid[x][y])/dx**2 + (Grid[x][idy0] + Grid[x][idyN] - 2*Grid[x][y])/dy**2;
    return res
}

void SAVE(int step, int Ngr)
{
    ofstream myfile("PhaseField_" + to_string(step) + "_" + to_string(Ngr+1)+ ".dat");
    for (int x=1; x<Nx+1; x++) {
        for (int y=1; y<Ny+1; y++) {
        myfile << double(ETAS[Ngr][x][y]) << "  ";
        }
        myfile << "\n";
    }
    myfile.close();
}
