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
// Physical constants
const double e  = 1.602177e-19; // [C]
const double kb = 1.380649e-23/e; // [eV/K]

// Discretization parameters
const int DIM = 3;         // N-dimension
const int MAXIT = 1000;    // Iteration number
const int NSTEP =  100;    // Step size for save
const int NMAX=  10;        // The number of grains
const int Nx = 100;        // Discretization x-direction
const int Ny = 100;        // Discretization y-direction
const int Nz =   1;        // Discretization y-direction
double dx, dy, dz = 1;     // Grid spacing in x-,y-,z-direction [Ang]
double dt = 0.005;         // Size of time step [ns]
double Mob=  10.0;         // Mobility [Angs^3/(eV/ns)]
double grad_coeff = 0.1;       // gradient coefficient
double A = 1, B = 1, Gam = 1;
const int BCells = 1; 	   // Number of boundary cells

// Phase-field storage
const int Mx = Nx + 2 * BCells;
const int My = Ny + 2 * BCells;
const int Mz = Nz + 2 * BCells;

double ETAS_TOT [Mx][My][Mz] = {0};
double ETAS     [NMAX][Mx][My][Mz] = {0};
double ETAS_new [NMAX][Mx][My][Mz] = {0};
double df         = 0;
double Lap_ETAS   = 0;
double SUM_eta_sq = 0;

// Functions
int check_stab(int dim);
double laplacian(double Grid[Mx][My][Mz], int x, int y, int z, double dx, double dy, double dz, int Nx, int Ny, int Nz);
double COUNT_AMOR(double Grid[Mx][My][Mz]);
double PROB_NUC(double T, double Va)
double growth_velocity(double T);
void SAVE(int step, int NMAX);

int main(int argc, char **argv)
{
    int Ngr = 0;
    int Ca;
    double Va, Xa;
    
    
    int rank, size, psize;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    psize = Nx / size;

    if (!check_stab(DIM)){
        MPI_Finalize();
        return 1;
    }

    // Initiallization
    srand(230917);
    int xc = rand() % (Nx) + 1;
    int yc = rand() % (Ny) + 1;
    int zc = rand() % (Nz) + 1;
    for (int x=1; x<Nx+1; x++) {
        for (int y=1; y<Ny+1; y++) {
            for (int z=1; z<Nz+1; z++) {
                if (sqrt((xc-x)*(xc-x) + (yc-y)*(yc-y) + (zc-z)*(zc-z)) < 1.2)
                ETAS[0][x][y][z] = 1.0;
            }
        } 
    }
    ETAS_TOT = ETAS_TOT + ETAS[0];
    
    // Evolve
    for (int k=0; k <= MAXIT; k++) {
        // Save data
        if (k % NSTEP == 0 && rank == 0){
            printf("%d/%d\n",k,MAXIT);
            for (int n=0; n < NMAX; n++){
                SAVE(k, n);
            }
        }
        
        Ca = COUNT_AMOR(ETAS_TOT);
        Va = Ca*dx*dy*dz;
        Xa = Ca/Nx*Ny*Nz;
        if (PROB_NUC(T, Va) > rand() && Xa < 0.9 && Ngr < NMAX)
        {
            Ngr++;
            xc, yc, zc = rand() % (Nx) + 1, rand() % (Ny) + 1, rand() % (Nz) + 1;
            while (ETAS_TOT[xc][yc][zc] > 0.1){
                xc, yc, zc = rand() % (Nx) + 1, rand() % (Ny) + 1, rand() % (Nz) + 1;
            }
            for (int x=1; x<Nx+1; x++) {
                for (int y=1; y<Ny+1; y++) {
                    for (int z=1; z<Nz+1; z++) {
                    if (sqrt((xc-x)*(xc-x) + (yc-y)*(yc-y) + (zc-z)*(zc-z)) < 1.2)
                    ETAS[Ngr][x][y][z] = 1.0;
                    }
                } 
            }
        }

        
    }
    /*
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
            if (rank == size -1){
                for (int x=rank*psize+1; x<Nx+1; x++)
                    MPI_Send(ETAS_new[N][x], sizeof(ETAS_new[N][x])/sizeof(double), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            } else if (rank != size -1 && rank != 0){
                for (int x=rank*psize+1; x<(rank+1)*psize+1; x++)
                    MPI_Send(ETAS_new[N][x], sizeof(ETAS_new[N][x])/sizeof(double), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            } else if (rank == 0){
                for (int r=1; r <  size; r++){
                    if (r == size -1){
                        for (int x=r*psize+1; x<Nx+1; x++)
                        MPI_Recv(ETAS_new[N][x], sizeof(ETAS_new[N][x])/sizeof(double), MPI_DOUBLE, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    } else {
                        for (int x=r*psize+1; x<(r+1)*psize+1; x++)
                        MPI_Recv(ETAS_new[N][x], sizeof(ETAS_new[N][x])/sizeof(double), MPI_DOUBLE, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }
            for (int x=1; x<Nx+1; x++)
            MPI_Bcast(ETAS_new[N][x], sizeof(ETAS_new[N][x])/sizeof(double), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
        } // NMAX loop


        // Periodic boundary condition
        for (int N=0; N < NMAX; N++){
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
        for (int N=0; N < NMAX; N++)
        for (int x=0; x<Mx; x++) 
        for (int y=0; y<My; y++)
            ETAS[N][x][y] = ETAS_new[N][x][y];
    } // MAXIT loop
    */
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    
    return 0; 
}

int check_stab(int dim)
{
    float res = (Mob*dt)*(1/(dx*dx)+1/(dy*dy)+1/(dz*dz));
    float crit= 1/2/dim
    if (res < crit) {
        printf("Maybe stable. %.2f is less than %.2f.\n", res, crit);
        return 1;
    } else {
        printf("Maybe unstable. %.2f is greater than %.2f.\n", res, crit);
        return 0;
    }
}

double laplacian(double Grid[Mx][My][Mz], int x, int y, int z, double dx, double dy, double dz, int Nx, int Ny, int Nz);
{
    int idx0, idxN = x - 1, x + 1;
    int idy0, idyN = y - 1, y + 1;
    int idz0, idzN = z - 1, z + 1;
    return (Grid[idx0][y][z] + Grid[idxN][y][z] - 2*Grid[x][y][z])/(dx*dx) + (Grid[x][idy0][z] + Grid[x][idyN][z] - 2*Grid[x][y][z])/(dy*dy) + (Grid[x][y][idz0] + Grid[x][y][idzN] - 2*Grid[x][y][z])/(dz*dz);
}

double COUNT_amor(double Grid[Mx][My][Mz])
{
    int count = 0;
    for (int x=1; x<Nx+1; x++) {
        for (int y=1; y<Ny+1; y++) {
            for (int z=1; z<Nz+1; z++) {
                if (Grid[x][y][z] < 0.5)
                {
                count++;
                }
            }
        }
    }
    return count
}

double PROB_NUC(double T, double Va)
{
    // Mater. Res. Soc. Symp. Proc. 989, 0989-A06-17 (2007)
    double beta   = 1/(kb*T);
    double I0 = 1.7e+44; // [m-3sec-1]
    double Ea = 5.3;     // [eV]
    double I_nuc = I0*exp(-beta*Ea);
    return 1-exp(-I_nuc*Va*dt);
}

double growth_velocity(double T)
{
    // Mater. Res. Soc. Symp. Proc. 989, 0989-A06-17 (2007)
    double beta   = 1/(kb*T);
    double u0 = 2.1e+7; // [m/s]
    double Ea = 3.1;    // [eV]
    return u0*exp(-beta*Ea);
}

void SAVE(int step, int NMAX)
{
    ofstream myfile("PhaseField_" + to_string(step) + "_" + to_string(NMAX)+ ".dat");
    for (int x=1; x<Nx+1; x++) {
        for (int y=1; y<Ny+1; y++) {
        myfile << fixed << setprecision(4) << double(ETAS[NMAX][x][y]);
        if (y < Ny) { myfile << "  ";}
        }
        myfile << "\n";
    }
    myfile.close();
}
