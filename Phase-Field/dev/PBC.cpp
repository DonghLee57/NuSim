// 'p': periodic, 'r': reflective, 'f': fixed
constexpr char PBCX = 'p'; 
constexpr char PBCY = 'p'; 
constexpr char PBCZ = 'p';

        // PBC setting
        cout << "Before" << endl;
        check_grid(ETAS_new[Ngr]);
        for (int N = 0; N < Ngr; ++N) {
            if (PBCX == 'p'){
                #pragma omp parallel for
                for (int y = 1; y < Ny + 1; ++y) {
                for (int z = 1; z < Nz + 1; ++z) {
                    ETAS_new[N][0][y][z] = ETAS_new[N][Nx][y][z];
                    ETAS_new[N][Mx][y][z] = ETAS_new[N][1][y][z];
                }}
            } else if (PBCX == 'r'){
                #pragma omp parallel for
                for (int y = 1; y < Ny + 1; ++y) {
                for (int z = 1; z < Nz + 1; ++z) {
                    ETAS_new[N][0][y][z] = ETAS_new[N][1][y][z];
                    ETAS_new[N][Mx][y][z] = ETAS_new[N][Nx][y][z];
                }}
            }
            if (PBCY == 'p'){
                #pragma omp parallel for
                for (int x = 1; x < Nx + 1; ++x) {
                for (int z = 1; z < Nz + 1; ++z) {
                    ETAS_new[N][x][0][z] = ETAS_new[N][x][Ny][z];
                    ETAS_new[N][x][My][z] = ETAS_new[N][x][1][z];
                }}
            } else if (PBCY == 'r'){
                #pragma omp parallel for
                for (int x = 1; x < Nx + 1; ++x) {
                for (int z = 1; z < Nz + 1; ++z) {
                    ETAS_new[N][x][0][z] = ETAS_new[N][x][1][z];
                    ETAS_new[N][x][My][z] = ETAS_new[N][x][Ny][z];
                }}
            }
            if (PBCZ == 'p'){
                #pragma omp parallel for
                for (int x = 1; x < Nx + 1; ++x) {
                for (int y = 1; y < Ny + 1; ++y) {
                    ETAS_new[N][x][y][0] = ETAS_new[N][x][y][Nz];
                    ETAS_new[N][x][y][Mz] = ETAS_new[N][x][y][1];
                }}
            } else if (PBCZ == 'r'){
                #pragma omp parallel for
                for (int x = 1; x < Nx + 1; ++x) {
                for (int y = 1; y < Ny + 1; ++y) {
                    ETAS_new[N][x][y][0] = ETAS_new[N][x][y][1];
                    ETAS_new[N][x][y][Mz] = ETAS_new[N][x][y][Nz];
            }}
            } 
        }
