const double PI = M_PI;
const double kb_eV = 8.617333262E-5; // [eV/K]
const double q = 1.602176634E-19;    // [C]
const double Tmelt = 1685;           // [K]
const double Vm = 1.205883199E-5;    // [m3/mol]
const double GAMMA = 0.38/q;         // Interfacial energy [eV/m2]
const double delta = 2E-10;          // Interfacial width  [m]
const double b = 2.2;                // ~2.2
const double dHf = 50500/q/Vm;       // [eV/m3]
const double K = 12*delta*GAMMA/b;   // Gradient coefficient [eV/m]
const double W = 6*GAMMA*b/delta;    // [eV/m3]

// Simulation parameters
const double Tsim = 1273;    // [K]
const int MAXIT =   5000;    // Iteration number
const int NSTEP =    500;    // Step size for save
const int Ngr=   1;          // The number of grains
const int Nx = 300;          // Discretization x-direction
const int Ny = 300;          // Discretization y-direction
double dx = 1E-10;           // Grid spacing in x-direction [m]
double dy = 1E-10;           // Grid spacing in y-direction [m]
double dt =   1E-7;           // Size of time step [s]
double Mob=  2.8E-24;          // Mobility [m^3/(eV/s)]
// ~8 nm/s
//double dt =   1E-4;           // Size of time step [s]
//double Mob=  1E-27;          // Mobility [m^3/(eV/s)]
double grad_coeff = K;       // gradient coefficient
double A = 1, B = 1;
const int BCells = 1;        // Number of boundary cells
