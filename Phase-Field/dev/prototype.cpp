#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <omp.h>
using namespace std;

constexpr double PI = M_PI;
constexpr double q = 1.602177e-19;         // [C]
constexpr double kb_eV = 1.380649e-23 / q; // [eV/K]
constexpr double Tmelt = 1685;             // [K]
constexpr double Vm = 1.205883199E-5;      // [m3/mol]
constexpr double GAMMA = 0.38/q;           // Interfacial energy [eV/m2]
constexpr double delta = 2E-10;            // Interfacial width  [m]
constexpr double b = 2.2;                  // ~2.2
constexpr double dHf = 50500/q/Vm;         // [eV/m3]
constexpr double K = 12*delta*GAMMA/b;     // Gradient coefficient [eV/m]
constexpr double W = 6*GAMMA*b/delta;      // [eV/m3]

const double Tsim = 1273; // [K]
const double Rc = 19; // [grid point]
const int DIM = 3;
const int MAXIT =    10000;
const int NSTEP =    1;
const int NMAX  = 4;
const int Nx = 8;
const int Ny = 8;
const int Nz = 1;
const int BCells = 1;

const double dx = 1E-10, dy = 1E-10, dz = 1E-10;  // Grid spacing in x-,y-,z-direction [m]
const double dt = 100E-9;                         // Size of time step [s]
const double Mob = 2.8E-24;                       // Mobility [m^3/(eV/s)]
double Lap_ETAS = 0;
double df= 0;

// Phase-field storage
const int Mx = Nx + 2 * BCells;
const int My = Ny + 2 * BCells;
const int Mz = Nz + 2 * BCells;

using PhaseField3D = vector<vector<vector<double>>>;

PhaseField3D ETAS_TOT(Mx, vector<vector<double>>(My, vector<double>(Mz, 0.0)));
vector<PhaseField3D> ETAS(NMAX, ETAS_TOT);
vector<PhaseField3D> ETAS_new(NMAX, ETAS_TOT);

// Function declarations
int CHECK_STABILITY(int dim);
double LAPLACIAN(const PhaseField3D& Grid, int x, int y, int z);
int COUNT_AMOR(const PhaseField3D& Grid);
double PROB_NUC(double T, double Va);
double GR_VEL(double T);
void SAVE_2D(int step, int NMAX, int z);
void check_grid(const PhaseField3D& Grid);
int SORT_GRAIN(vector<PhaseField3D>& NGrid, int Ngr);

int main()
{
    int Ngr = 0;
    int Ca;
    double Va, Xa;
    double distance;
    double SUM_eta_sq = 0;

    //if (!CHECK_STABILITY(DIM)) {
    //    return 1;
    //}

    // Initialization
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distrX(1, Nx);
    uniform_int_distribution<int> distrY(1, Ny);
    uniform_int_distribution<int> distrZ(1, Nz);
    uniform_real_distribution<double> probDist(0.0, 1.0);

    int rx = distrX(gen);
    int ry = distrY(gen);
    int rz = distrZ(gen);
    int z0 = rz;
    #pragma omp parallel for private(distance)
    for (int x = 1; x < Nx + 1; ++x) {
        for (int y = 1; y < Ny + 1; ++y) {
            for (int z = 1; z < Nz + 1; ++z) {
                distance = sqrt((rx - x) * (rx - x) + (ry - y) * (ry - y) + (rz - z) * (rz - z));
                ETAS[Ngr][x][y][z] = (distance < Rc) ? 1.0 : 0.0;
            }
        }
    }
    ETAS_TOT = ETAS[Ngr];

    Ngr++;
    // Evolve
    for (int k = 0; k <= MAXIT; ++k) {
        
        // Save data
        if (k % NSTEP == 0) {
            cout << k << '/' << MAXIT << '\n';
            for (int n = 0; n < NMAX; ++n) {
                SAVE_2D(k, n, z0);
            }
        }
	
        // Test
        Ca = COUNT_AMOR(ETAS_TOT);
        Va = Ca * dx * dy * dz;
        Xa = (double) Ca / (Nx * Ny * Nz);
        cout << "No. of Amor.: " << Ca << endl;
        cout << "Amor. Vol.: " << Va << " (m^3)"  << endl;
        cout << "Amor. Fraction: " << Xa << endl;

        if (PROB_NUC(Tsim, Va) > probDist(gen) && Xa > 0.1 && Ngr < NMAX) {
            rx = distrX(gen);
            ry = distrY(gen);
            rz = distrZ(gen);
            while (ETAS_TOT[rx][ry][rz] > 0.1) {
                rx = distrX(gen);
                ry = distrY(gen);
                rz = distrZ(gen);
            }
            #pragma omp parallel for private(distance)
            for (int x = 1; x < Nx + 1; ++x) {
                for (int y = 1; y < Ny + 1; ++y) {
                    for (int z = 1; z < Nz + 1; ++z) {
                        distance = sqrt((rx - x) * (rx - x) + (ry - y) * (ry - y) + (rz - z) * (rz - z));
                        ETAS[Ngr][x][y][z] = (distance < Rc) ? 1.0 : 0.0;
			ETAS_TOT[x][y][z] = ETAS_TOT[x][y][z] + ETAS[Ngr][x][y][z];
                    }
                }
            }
            Ngr++;
        }

        // Calculate
        for (int N = 0; N < Ngr; ++N) {
        #pragma omp parallel for private(SUM_eta_sq, df, Lap_ETAS)
        for (int x = 1; x < Nx + 1; ++x) {
            for (int y = 1; y < Ny + 1; ++y) {
                for (int z = 1; z < Nz + 1; ++z) {
                    SUM_eta_sq = 0;
                    for (int ngr = 0; ngr < Ngr; ++ngr) {
                        SUM_eta_sq += pow(ETAS[ngr][x][y][z], 2);
                    }
                    df = 30 * pow(ETAS[N][x][y][z],2) * pow(1-ETAS[N][x][y][z],2) * (-dHf * (1 - Tsim/Tmelt))
                       +  2 * W * ETAS[N][x][y][z] * (1 - ETAS[N][x][y][z]) * (1 - 2*ETAS[N][x][y][z])
                       + 2E+29 * 2 * ETAS[N][x][y][z] * (SUM_eta_sq - pow(ETAS[N][x][y][z],2));
                    Lap_ETAS = LAPLACIAN(ETAS[N], x, y, z);
                    ETAS_new[N][x][y][z] = ETAS[N][x][y][z] - dt * Mob * (df - K * Lap_ETAS);
                }
            }
        }
        }

        // Update ETAS
        //ETAS = ETAS_new;
        for (int N = 0; N < Ngr; ++N) {
        #pragma omp parallel for
        for (int x = 1; x < Nx + 1; ++x) {
        for (int y = 1; y < Ny + 1; ++y) {
        for (int z = 1; z < Nz + 1; ++z) {
            if (N == 0){
                ETAS_TOT[x][y][z] = ETAS_new[N][x][y][z];
            } else {
                ETAS_TOT[x][y][z] += ETAS_new[N][x][y][z];
            }
            ETAS[N][x][y][z] = ETAS_new[N][x][y][z];
        }}}}
        // position
        Ngr = SORT_GRAIN_GRID(ETAS, Ngr);
    } // MAXIT loop
    return 0;
}

int CHECK_STABILITY(int dim)
{
    float res = (Mob * dt) * (1 / (dx * dx) + 1 / (dy * dy) + 1 / (dz * dz));
    float crit = 1.0 / 2.0 / dim;
    if (res < crit) {
        cout << "Maybe stable. " << res << " is less than " << crit << ".\n";
        return 1;
    } else {
        cout << "Maybe unstable. " << res << " is greater than " << crit << ".\n";
        return 0;
    }
}

double LAPLACIAN(const PhaseField3D& Grid, int x, int y, int z)
{
    int idx0 = x - 1, idxN = x + 1;
    int idy0 = y - 1, idyN = y + 1;
    int idz0 = z - 1, idzN = z + 1;

    return (Grid[idx0][y][z] + Grid[idxN][y][z] - 2 * Grid[x][y][z]) / (dx * dx) +
           (Grid[x][idy0][z] + Grid[x][idyN][z] - 2 * Grid[x][y][z]) / (dy * dy) +
           (Grid[x][y][idz0] + Grid[x][y][idzN] - 2 * Grid[x][y][z]) / (dz * dz);
}

int COUNT_AMOR(const PhaseField3D& Grid)
{
    int count = 0;
    #pragma omp parallel for reduction(+ : count)
    for (int x = 1; x < Nx + 1; ++x) {
        for (int y = 1; y < Ny + 1; ++y) {
            for (int z = 1; z < Nz + 1; ++z) {
                if (Grid[x][y][z] < 0.5) {
                    ++count;
                }
            }
        }
    }
    return count;
}

double PROB_NUC(double T, double Va)
{
    // Mater. Res. Soc. Symp. Proc. 989, 0989-A06-17 (2007)
    double beta = 1 / (kb_eV * T);
    double I0 = 1.7e+44; // [m-3sec-1]
    double Ea = 5.3;     // [eV]
    double I_nuc = I0 * exp(-beta * Ea);
    return 1 - exp(-I_nuc * Va * dt);
}

double GR_VEL(double T)
{
    // Mater. Res. Soc. Symp. Proc. 989, 0989-A06-17 (2007)
    double beta = 1 / (kb_eV * T);
    double u0 = 2.1e+7; // [m/s]
    double Ea = 3.1;    // [eV]
    return u0 * exp(-beta * Ea);
}

void SAVE_2D(int step, int NMAX, int z)
{
    ofstream myfile("PhaseField_" + to_string(step) + "_" + to_string(NMAX) + ".dat");
    for (int x = 1; x < Nx + 1; ++x) {
        for (int y = 1; y < Ny + 1; ++y) {
            myfile << fixed << setprecision(4) << ETAS[NMAX][x][y][z];
            if (y < Ny) {
                myfile << "  ";
            }
        }
        myfile << "\n";
    }
    myfile.close();
}

int SORT_GRAIN(vector<PhaseField3D>& NGrid, int Ngr){
    vector<PhaseField3D> sorted(NMAX, ETAS_TOT);
    int sorted_Ngr = 0;
    int n_amor;
    int crit = 100;
    for (int N = 0; N < Ngr; N++){
        n_amor = COUNT_AMOR(NGrid[N]);
        if ( (Nx*Ny*Nz - n_amor) > crit){
            sorted[sorted_Ngr] = NGrid[N];
            sorted_Ngr++;
        }
    }
    ETAS = sorted;
    return sorted_Ngr;
}
