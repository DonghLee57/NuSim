#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <random>
#include <mpi.h>
using namespace std;

// Discretization parameters
const int MAXIT = 1000;    // Iteration number
const int NSTEP =  100;    // Step size for save
const int Ngr=  10;        // The number of grains
const int Nx = 100;        // Discretization x-direction
const int Ny = 100;        // Discretization y-direction
double dx = 1;             // Grid spacing in x-direction [Ang]
double dy = 1;             // Grid spacing in y-direction [Ang]
double dt = 0.005;          // Size of time step [ns]
double Mob=  10.0;             // Mobility [Angs^3/(eV/ns)]
double grad_coeff = 0.1;       // gradient coefficient
double A = 1, B = 1, Gam = 1;
const int BCells = 1; 	       // Number of boundary cells

// Phase-field storage
const int Mx = Nx + 2 * BCells;
const int My = Ny + 2 * BCells;

double ETAS     [Ngr][Mx][My] = {0};
double ETAS_new [Ngr][Mx][My] = {0};
double df         = 0;
double Lap_ETAS   = 0;
double SUM_eta_sq = 0;

// Functions
int check_stab();
double lap_2D(double Grid[Mx][My], int x, int y, double dx, double dy, int Nx, int Ny);
void SAVE(int step, int Ngr);

int main(int argc, char **argv)
{
    int rank, size, psize;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    psize = Nx / size;

    if (!check_stab()){
        MPI_Finalize();
        return 1;
    }

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

    // Evolve
    for (int k=0; k <= MAXIT; k++) {
        // Save data
        if (k % NSTEP == 0 && rank == 0){
            printf("%d/%d\n",k,MAXIT);
            for (int n=0; n < Ngr; n++){
                SAVE(k, n);
            }
        }
        // Calculate
        for (int N=0; N < Ngr; N++){
            if (rank == size - 1){
                for (int x=rank*psize+1; x<Nx+1; x++) {
                for (int y=1; y<Ny+1; y++) {
                    SUM_eta_sq = 0;
                    for (int ngr=0; ngr < Ngr; ngr++) SUM_eta_sq = SUM_eta_sq + pow(ETAS[ngr][x][y],2);
                    df = -A*ETAS[N][x][y] + B*pow(ETAS[N][x][y],3) + 2*ETAS[N][x][y]*(SUM_eta_sq - pow(ETAS[N][x][y],2));
                    Lap_ETAS = lap_2D(ETAS[N], x, y, dx, dy, Nx, Ny);
                    ETAS_new[N][x][y] = ETAS[N][x][y] - dt*Mob*(df - grad_coeff*Lap_ETAS);
                }
                }
            } else {
                for (int x=rank*psize+1; x<(rank+1)*psize+1; x++) {
                for (int y=1; y<Ny+1; y++) {
                    SUM_eta_sq = 0;
                    for (int ngr=0; ngr < Ngr; ngr++) SUM_eta_sq = SUM_eta_sq + pow(ETAS[ngr][x][y],2);
                    df = -A*ETAS[N][x][y] + B*pow(ETAS[N][x][y],3) + 2*ETAS[N][x][y]*(SUM_eta_sq - pow(ETAS[N][x][y],2));
                    Lap_ETAS = lap_2D(ETAS[N], x, y, dx, dy, Nx, Ny);
                    ETAS_new[N][x][y] = ETAS[N][x][y] - dt*Mob*(df - grad_coeff*Lap_ETAS);
                }
                }
            }
            // MPI comm
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == size-1){
                for (int x=rank*psize+1; x<Nx+1; x++) 
                MPI_Send(ETAS_new[N][x], sizeof(ETAS_new[N][x])/sizeof(double), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            } else {
                for (int x=rank*psize+1; x<(rank+1)*psize+1; x++) 
                MPI_Send(ETAS_new[N][x], sizeof(ETAS_new[N][x])/sizeof(double), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }

            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == 0) {
                for (int r = 1; r < size; r++) {
                    if (r == size-1) {
                        for (int x=r*size+1; x<Nx+1; x++)
                        MPI_Recv(ETAS_new[N][x], sizeof(ETAS_new[N][x])/sizeof(double), MPI_DOUBLE, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    } else {
                        for (int x=r*size+1; x<(r+1)*size+1; x++)
                        MPI_Recv(ETAS_new[N][x], sizeof(ETAS_new[N][x])/sizeof(double), MPI_DOUBLE, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
                for (int n=0; n<Nx+1; n++)
                MPI_Bcast(ETAS_new[n], sizeof(ETAS_new[n])/sizeof(double), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);
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
    } // MAXIT loop

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

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
    for (int x=1; x<Nx+1; x++) {
        for (int y=1; y<Ny+1; y++) {
        myfile << fixed << setprecision(4) << double(ETAS[Ngr][x][y]);
        if (y < Ny) { myfile << "  ";}
        }
        myfile << "\n";
    }
    myfile.close();
}

/*
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <random>
#include <mpi.h>
using namespace std;

int main(int argc, char **argv)
{
    int rank, size, psize;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    psize = 4 / size;

    int MAT[4][4] = {0};

    for (int n=rank*psize; n<(rank+1)*psize; n++){
    printf(" %d ", rank);
    for (int m=0; m<4; m++){
        MAT[n][m] = rank;
        printf(" %d ", MAT[n][m]);
    }
    printf(" \n ");
    }

    MPI_Barrier(MPI_COMM_WORLD);
    printf(" \n ");
    
    //for (int r = 1; r < size; r++) {
        //for (int n=r*size; n<(r+1)*size; n++){
        for (int n=rank*size; n<(rank+1)*size; n++){
            MPI_Send(MAT[n], sizeof(MAT[n])/sizeof(int), MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    //}
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank ==0){
    for (int r = 1; r < size; r++) {
        for (int n=r*size; n<(r+1)*size; n++){
        MPI_Recv(MAT[n], sizeof(MAT[n])/sizeof(int), MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    for (int n=0; n<4; n++){
        MPI_Bcast(MAT[n], sizeof(MAT[n])/sizeof(int), MPI_INT, 0, MPI_COMM_WORLD);
    }
    }
    MPI_Barrier(MPI_COMM_WORLD);


    for (int n=0; n<4; n++){
    printf(" %d ", rank);
    for (int m=0; m<4; m++){
        printf(" %d ", MAT[n][m]);
    }
    printf(" \n ");
    }


    MPI_Finalize();
    return 0; 
}
*/
