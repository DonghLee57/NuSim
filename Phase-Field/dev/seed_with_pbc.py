            xc = np.random.randint(low=0,high=Nx)
            yc = np.random.randint(low=0,high=Ny)
            for x in range(Nx):
                for y in range(Ny):
                    if x0[0] == 'p':
                        if np.fabs(xc-x) > Nx//2 and np.sign(xc-x) < 0: x -= Nx
                    elif x0[0] == 'f':
                        if np.fabs(xc-x) > Nx//2 and np.sign(xc-x) < 0: x = xc
                    if xN[0] == 'p':
                        if np.fabs(xc-x) > Nx//2 and np.sign(xc-x) > 0: x += Nx
                    elif xN[0] == 'f':
                        if np.fabs(xc-x) > Nx//2 and np.sign(xc-x) > 0: x = xc
                    if y0[0] == 'p':
                        if np.fabs(yc-y) > Ny//2 and np.sign(yc-y) < 0: y -= Ny
                    elif y0[0] == 'f':
                        if np.fabs(yc-y) > Nx//2 and np.sign(yc-y) < 0: y = yc
                    if yN[0] == 'p':
                        if np.fabs(yc-y) > Ny//2 and np.sign(yc-y) > 0: y += Ny
                    elif y0[0] == 'f':
                        if np.fabs(yc-y) > Nx//2 and np.sign(yc-y) > 0: y = yc
                    if ((xc - x)**2 + (yc - y)**2)**0.5 < R:
                        ETAS[N][x%Nx][y%Ny] = 1.0
