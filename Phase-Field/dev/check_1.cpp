#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <random>
using namespace std;
// Physical constants
constexpr double e     = 1.602177e-19;   // [C]
constexpr double kb_eV = 1.380649e-23/e; // [eV/K]

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

double ETAS_TOT [Mx][My][Mz];
double ETAS     [NMAX][Mx][My][Mz];
double ETAS_new [NMAX][Mx][My][Mz];
double df         = 0;
double Lap_ETAS   = 0;
double SUM_eta_sq = 0;

// Functions
int CHECK_STABILITY(int dim);
double LAPLACIAN(double Grid[Mx][My][Mz], int x, int y, int z, double dx, double dy, double dz, int Nx, int Ny, int Nz);
double COUNT_AMOR(double Grid[Mx][My][Mz]);
double PROB_NUC(double T, double Va);
double GR_VEL(double T);
void SAVE_2D(int step, int NMAX);

int main(int argc, char **argv)
{
    int Ngr = 0;
    int Ca;
    double Va, Xa;

    if (!CHECK_STABILITY(DIM)){
        return 1;
    }

    // Initiallization
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distrX(1,Nx);
    uniform_int_distribution<int> distrY(1,Ny);
    uniform_int_distribution<int> distrZ(1,Nz);
    int rx = distrX(gen);
    int ry = distrY(gen);
    int rz = distrZ(gen);
    for (int x=1; x<Nx+1; x++) {
        for (int y=1; y<Ny+1; y++) {
            for (int z=1; z<Nz+1; z++) {
                if (sqrt((rx-x)*(rx-x) + (ry-y)*(ry-y) + (rz-z)*(rz-z)) < 1.2)
                    ETAS[0][x][y][z] = 1.0;
                else
                    ETAS[0][x][y][z] = 0.0;
            }
        } 
    }
    
    copy(&ETAS[0][0][0][0], &ETAS[0][0][0][0] + Mx*My*Mz, &ETAS_TOT[0][0][0]);
    
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
        Xa = Ca/(Nx*Ny*Nz);
        if (PROB_NUC(T, Va) > uniform_real_distribution<>(0.0, 1.0)(gen) && Xa < 0.9 && Ngr < NMAX)
        {
            Ngr++;
            rx = distrX(gen);
            ry = distrY(gen);
            rz = distrZ(gen);
            while (ETAS_TOT[rx][ry][rz] > 0.1){
                rx = distrX(gen);
                ry = distrY(gen);
                rz = distrZ(gen);
            }
            for (int x=1; x<Nx+1; x++) {
                for (int y=1; y<Ny+1; y++) {
                    for (int z=1; z<Nz+1; z++) {
                    if (sqrt((rx-x)*(rx-x) + (ry-y)*(ry-y) + (rz-z)*(rz-z)) < 1.2)
                    ETAS[Ngr][x][y][z] = 1.0;
                    }
                } 
            }
        }

        // Calculate
        for (int N=0; N < Ngr; N++){
            for (int x=rank*psize+1; x<Nx+1; x++) {
                for (int y=1; y<Ny+1; y++) {
                    SUM_eta_sq = 0;
                    for (int ngr=0; ngr < Ngr; ngr++) SUM_eta_sq = SUM_eta_sq + pow(ETAS[ngr][x][y][z],2);
                    df = -A*ETAS[N][x][y][z] + B*pow(ETAS[N][x][y][z],3) + 2*ETAS[N][x][y][z]*(SUM_eta_sq - pow(ETAS[N][x][y][z],2));
                    Lap_ETAS = LAPLACIAN(ETAS[N], x, y, z);
                    ETAS_new[N][x][y][z] = ETAS[N][x][y][z] - dt*Mob*(df - grad_coeff*Lap_ETAS);
                }
            }

        }

        // ETAS 
        for (int N=0; N < NMAX; N++)
        for (int x=0; x<Mx; x++) 
        for (int y=0; y<My; y++)
        for (int z=0; z<Mz; z++)
            ETAS[N][x][y][z] = ETAS_new[N][x][y][z];
    } // MAXIT loop
    
    return 0; 
}

int CHECK_STABILITY(int dim)
{
    float res = (Mob*dt)*(1/(dx*dx)+1/(dy*dy)+1/(dz*dz));
    float crit= 1/2/dim;
    if (res < crit) {
        printf("Maybe stable. %.2f is less than %.2f.\n", res, crit);
        return 1;
    } else {
        printf("Maybe unstable. %.2f is greater than %.2f.\n", res, crit);
        return 0;
    }
}

double LAPLACIAN(double Grid[Mx][My][Mz], int x, int y, int z);
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
    double beta = 1/(kb_eV*T);
    double I0 = 1.7e+44; // [m-3sec-1]
    double Ea = 5.3;     // [eV]
    double I_nuc = I0*exp(-beta*Ea);
    return 1-exp(-I_nuc*Va*dt);
}

double GR_VEL(double T)
{
    // Mater. Res. Soc. Symp. Proc. 989, 0989-A06-17 (2007)
    double beta   = 1/(kb_eV*T);
    double u0 = 2.1e+7; // [m/s]
    double Ea = 3.1;    // [eV]
    return u0*exp(-beta*Ea);
}

void SAVE_2D(int step, int NMAX, int z)
{
    ofstream myfile("PhaseField_" + to_string(step) + "_" + to_string(NMAX)+ "_" + to_string(z)+ ".dat");
    for (int x=1; x<Nx+1; x++) {
        for (int y=1; y<Ny+1; y++) {
        myfile << fixed << setprecision(4) << double(ETAS[NMAX][x][y][z]);
        if (y < Ny) { myfile << "  ";}
        }
        myfile << "\n";
    }
    myfile.close();
}
