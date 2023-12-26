#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>

constexpr double e = 1.602177e-19;   // [C]
constexpr double kb_eV = 1.380649e-23 / e; // [eV/K]

const int DIM = 3;
const int MAXIT = 1000;
const int NSTEP = 100;
const int NMAX = 10;
const int Nx = 100;
const int Ny = 100;
const int Nz = 1;
const int BCells = 1;

double dx = 1.0, dy = 1.0, dz = 1.0;   // Grid spacing in x-,y-,z-direction [Ang]
double dt = 0.001;                     // Size of time step [ns]
double Mob = 1.0;                     // Mobility [Angs^3/(eV/ns)]
double grad_coeff = 0.1;               // Gradient coefficient
double A = 1.0, B = 1.0, Gam = 1.0;
double Lap_ETAS = 0;
double df= 0;

// Phase-field storage
const int Mx = Nx + 2 * BCells;
const int My = Ny + 2 * BCells;
const int Mz = Nz + 2 * BCells;

using PhaseField3D = std::vector<std::vector<std::vector<double>>>;

PhaseField3D ETAS_TOT(Mx, std::vector<std::vector<double>>(My, std::vector<double>(Mz, 0.0)));
std::vector<PhaseField3D> ETAS(NMAX, ETAS_TOT);
std::vector<PhaseField3D> ETAS_new(NMAX, ETAS_TOT);

// Function declarations
int CHECK_STABILITY(int dim);
double LAPLACIAN(const PhaseField3D& Grid, int x, int y, int z);
int COUNT_AMOR(const PhaseField3D& Grid);
double PROB_NUC(double T, double Va);
double GR_VEL(double T);
void SAVE_2D(int step, int NMAX, int z);

int main()
{
    int Ngr = 0;
    int Ca;
    double Va, Xa;

    if (!CHECK_STABILITY(DIM)) {
        return 1;
    }

    // Initialization
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distrX(1, Nx);
    std::uniform_int_distribution<int> distrY(1, Ny);
    std::uniform_int_distribution<int> distrZ(1, Nz);

    int rx = distrX(gen);
    int ry = distrY(gen);
    int rz = distrZ(gen);
    for (int x = 1; x < Nx + 1; ++x) {
        for (int y = 1; y < Ny + 1; ++y) {
            for (int z = 1; z < Nz + 1; ++z) {
                double distance = std::sqrt((rx - x) * (rx - x) + (ry - y) * (ry - y) + (rz - z) * (rz - z));
                ETAS[0][x][y][z] = (distance < 1.2) ? 1.0 : 0.0;
            }
        }
    }
    Ngr++;

    ETAS_TOT = ETAS[0];

    // Evolve
    for (int k = 0; k <= MAXIT; ++k) {
        // Save data
        if (k % NSTEP == 0) {
            std::cout << k << '/' << MAXIT << '\n';
            for (int n = 0; n < NMAX; ++n) {
                SAVE_2D(k, n, 0);
            }
        }

        Ca = COUNT_AMOR(ETAS_TOT);
        Va = Ca * dx * dy * dz;
        Xa = Ca / (Nx * Ny * Nz);

        std::uniform_real_distribution<double> probDist(0.0, 1.0);

        if (PROB_NUC(dt, Va) > probDist(gen) && Xa < 0.9 && Ngr < NMAX) {
            ++Ngr;
            rx = distrX(gen);
            ry = distrY(gen);
            rz = distrZ(gen);

            while (ETAS_TOT[rx][ry][rz] > 0.1) {
                rx = distrX(gen);
                ry = distrY(gen);
                rz = distrZ(gen);
            }

            for (int N=0; N < Ngr; N++) {
            for (int x = 1; x < Nx + 1; ++x) {
                for (int y = 1; y < Ny + 1; ++y) {
                    for (int z = 1; z < Nz + 1; ++z) {
                        double distance = std::sqrt((rx - x) * (rx - x) + (ry - y) * (ry - y) + (rz - z) * (rz - z));
                        ETAS[N][x][y][z] = (distance < 1.2) ? 1.0 : 0.0;
                    }
                }
            }
            }
        }

        // Calculate
        for (int N = 0; N < Ngr; ++N) {
        for (int x = 1; x < Nx + 1; ++x) {
            for (int y = 1; y < Ny + 1; ++y) {
                for (int z = 1; z < Nz + 1; ++z) {
                    double SUM_eta_sq = 0;
                    for (int ngr = 0; ngr < Ngr; ++ngr) {
                        SUM_eta_sq += std::pow(ETAS[ngr][x][y][z], 2);
                    }
                    df = -A * ETAS[N][x][y][z] + B * std::pow(ETAS[N][x][y][z], 3) +
                         2 * ETAS[N][x][y][z] * (SUM_eta_sq - std::pow(ETAS[N][x][y][z], 2));
                    Lap_ETAS = LAPLACIAN(ETAS[N], x, y, z);
                    ETAS_new[N][x][y][z] = ETAS[N][x][y][z] - dt * Mob * (df - grad_coeff * Lap_ETAS);
                }
            }
        }
        }
        // Update ETAS
        ETAS = ETAS_new;
    } // MAXIT loop

    return 0;
}

int CHECK_STABILITY(int dim)
{
    float res = (Mob * dt) * (1 / (dx * dx) + 1 / (dy * dy) + 1 / (dz * dz));
    float crit = 1.0 / 2.0 / dim;
    if (res < crit) {
        std::cout << "Maybe stable. " << res << " is less than " << crit << ".\n";
        return 1;
    } else {
        std::cout << "Maybe unstable. " << res << " is greater than " << crit << ".\n";
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
    std::ofstream myfile("PhaseField_" + std::to_string(step) + "_" + std::to_string(NMAX) + "_" + std::to_string(z) + ".dat");
    for (int x = 1; x < Nx + 1; ++x) {
        for (int y = 1; y < Ny + 1; ++y) {
            myfile << std::fixed << std::setprecision(4) << ETAS[NMAX][x][y][z];
            if (y < Ny) {
                myfile << "  ";
            }
        }
        myfile << "\n";
    }
    myfile.close();
}
