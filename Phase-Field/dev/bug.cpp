
int COUNT_AMOR(const PhaseField3D& Grid)
{
    int count = 0;
    #pragma omp parallel for reduction(+ : count)
    for (int x = 1; x < Nx + 1; ++x) {
        for (int y = 1; y < Ny + 1; ++y) {
            for (int z = 1; z < Nz + 1; ++z) {
                if (Grid[x][y][z] < 0.5) {
                    count = count + 1;
                }
            }
        }
    }
    return count;
}
