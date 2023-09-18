#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <random>
#include <omp.h>
using namespace std;

// Discretization parameters
const int MAXIT =  300;    // Iteration number
const int FROMIT = 0;      // Continue from FROMIT
const int NSTEP =  100;    // Step size for save
const int Ngr=  10;        // The number of grains
const int Nx = 100;        // Discretization x-direction
const int Ny = 100;        // Discretization y-direction
double dx = 1;             // Grid spacing in x-direction [Ang]
double dy = 1;             // Grid spacing in y-direction [Ang]
double dt = 0.01;          // Size of time step [s]
double Mob=  10.0;             // Mobility
double grad_coeff = 0.1;       // gradient coefficient
double A = 1, B = 1, Gam = 1;
const int BCells = 1; 	       // Number of boundary cells

// Phase-field storage
const int Mx = Nx + 2 * BCells;
const int My = Ny + 2 * BCells;

double ETAS     [Ngr][Mx][My] = {0};
double ETAS_new [Ngr][Mx][My] = {0};
double df       [Mx][My] = {0};
double Lap_ETAS   = 0;
double SUM_eta_sq = 0;
int    MAX_Ngr    = Ngr;

// Functions
int check_stab();
double lap_2D(double Grid[Mx][My], int x, int y, double dx, double dy, int Nx, int Ny);
void SAVE(int step, int Ngr);
void LOAD(int step, int Ngr);

int main()
{
    check_stab();

    // Initiallization
    srand(230917);
    for (int n=0; n < Ngr; n++) {
        int xc = rand() % (Nx) + 1;
        int yc = rand() % (Ny) + 1;
        for (int x=1; x<Nx+1; x++) {
            for (int y=1; y<Ny+1; y++) {
                if (sqrt(pow(xc-x,2) + pow(yc-y,2)) < 1.2)
                {
                    ETAS[n][x][y] = 1.0;
                }
            }
        } 
    }

    // Load data
    if (FROMIT > 0) {
        for (int n=0; n < Ngr; n++)
            LOAD(FROMIT, n);
    }

    // Evolve
    #pragma omp parallel
    //int N_THR = omp_get_num_threads();
    cout << to_string(omp_get_num_threads()) << " threads\n";
    for (int k=FROMIT+1; k <= FROMIT+MAXIT; k++) {
        // Save data
        if (k % NSTEP == 0){
            printf("%d/%d\n",k,FROMIT+MAXIT);
            for (int n=0; n < Ngr; n++){
            SAVE(k, n);
            }
        }
        // Calculate
        for (int N=0; N < Ngr; N++){
            for (int x=1; x<Nx+1; x++) {
            for (int y=1; y<Ny+1; y++) {
                SUM_eta_sq = 0;
                for (int ngr=0; ngr < Ngr; ngr++){
                    if (N != ngr) {SUM_eta_sq = SUM_eta_sq + pow(ETAS[ngr][x][y],2);}
                }
                df[x][y] = -A*ETAS[N][x][y] + B*pow(ETAS[N][x][y],3) + 2*ETAS[N][x][y]*SUM_eta_sq;
            }
            }
            for (int x=1; x<Nx+1; x++){
            for (int y=1; y<Ny+1; y++){ 
                Lap_ETAS = lap_2D(ETAS[N], x, y, dx, dy, Nx, Ny);
                ETAS_new[N][x][y] = ETAS[N][x][y] - dt*Mob*(df[x][y] - grad_coeff*Lap_ETAS);
            }
            }
        } // Ngr loop

        // Periodic boundary condition
        for (int N=0; N < Ngr; N++){
        for (int y=0; y < My; y++){
            ETAS_new[N][0][y]    =  ETAS_new[N][Nx][y];
            ETAS_new[N][Nx+1][y] =  ETAS_new[N][1][y];
        }
        for (int x=0; x < Mx; x++){
            ETAS_new[N][x][0]    =  ETAS_new[N][x][Ny];
            ETAS_new[N][x][Ny+1] =  ETAS_new[N][x][1];
        }
        }

        // ETAS 
        for (int N=0; N < Ngr; N++)
        for (int x=0; x<Mx; x++) 
        for (int y=0; y<My; y++)
            ETAS[N][x][y] = ETAS_new[N][x][y];
        #pragma omp barrier
    } // MAXIT loop

    return 0; 
}

int check_stab()
{
    float res = (Mob*dt)*(1/pow(dx,2)+1/pow(dy,2));
    if (res < 0.25) {
        printf("Maybe stable. %.2f is less than 0.25.\n", res);
        return 1;
    } else {
        printf("Maybe unstable. %.2f is greater than 0.25.\n", res);
        return 0;
    }
}

double lap_2D(double Grid[Mx][My], int x, int y, double dx, double dy, int Nx, int Ny)
{
    int idx0 = x - 1;
    int idxN = x + 1;
    int idy0 = y - 1;
    int idyN = y + 1;
    return (Grid[idx0][y] + Grid[idxN][y] - 2*Grid[x][y])/pow(dx,2) + (Grid[x][idy0] + Grid[x][idyN] - 2*Grid[x][y])/pow(dy,2);
}

void SAVE(int step, int Ngr)
{
    ofstream myfile("PhaseField_" + to_string(step) + "_" + to_string(Ngr)+ ".dat");
    int NSITE = 0;
    for (int x=1; x<Nx+1; x++) {
        for (int y=1; y<Ny+1; y++) {
        if (ETAS[Ngr][x][y] > 0.9) {NSITE = NSITE + 1;}
        myfile << double(ETAS[Ngr][x][y]);
        if (y < Ny) { myfile << "  ";}
        }
        myfile << "\n";
    }
    if (NSITE == 0) {
        MAX_Ngr = MAX_Ngr - 1;
        printf("Effective %d grains are in my system.\n", MAX_Ngr);
    }
    myfile.close();
}

void LOAD(int step, int Ngr)
{
    ifstream myfile("PhaseField_" + to_string(step) + "_" + to_string(Ngr)+ ".dat");
    int NSITE = 0;
    for (int x=1; x<Nx+1; x++) {
        for (int y=1; y<Ny+1; y++) {
        myfile >> ETAS[Ngr][x][y];
        if (ETAS[Ngr][x][y] > 0.9) {NSITE = NSITE + 1;}
        }
        myfile << "\n";
    }
    if (NSITE == 0) {
        MAX_Ngr = MAX_Ngr - 1;
        printf("Effective %d grains are in my system.\n", MAX_Ngr);
    }
    myfile.close();
}
