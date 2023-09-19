import numpy as np
import matplotlib.pyplot as plt

MAXIT = 1001
Nx = 100
Mx = Nx + 2
dx = 1
dt = 0.01
L = 1.0
kappa = 1.0
W = 1.0/4.0

def g(eta,W=W):
    return W*(eta)**2*(1-eta)**2

def deriv_g(eta,W=W):
    return 2*W*eta*(1-eta)*(1-2*eta)

def LAP_1D(eta, idx):
    return (eta[idx+1] + eta[idx-1] - 2*eta[idx])/dx**2

def analytic_sol(X, Xc, W=W, kappa=kappa):
    return 0.5*(1 + np.tanh(np.sqrt(W/(2*kappa))*(X-Xc)))

def VN_stability():
    eval = dt*L/dx**2
    print(f"VN_stability: {eval:.4f}")
    return eval

def evolve(L=L, kappa=kappa, W=W):
    ETA = np.zeros((MAXIT,Mx))
    for i in range(MAXIT-1):
        if i == 0:
             ETA[i][Mx//2:].fill(1)
             ETA[:,-1].fill(1)
        for x in range(1,Mx-1):
             ETA[i+1][x] = ETA[i][x] - dt*L*(deriv_g(ETA[i][x],W) - kappa*LAP_1D(ETA[i],x))
    return ETA

if VN_stability() < 0.5:
    # Analytic solution vs Finite difference
    fig, ax = plt.subplots(1, 1, figsize=(6,6))
    X = np.linspace(0,Mx,1000)
    ax.plot(evolve()[-1], label='Finite difference')
    ax.plot(X, analytic_sol(X,Mx//2), ls='--',label='Analytic solution')
    ax.set_xlabel('x',fontsize=12)
    ax.set_ylabel(r'$\eta$',fontsize=12)
    plt.legend()
    plt.savefig('FD_vs_Analytic.png')
    plt.close()

    # W dependence
    fig, ax = plt.subplots(1, 1, figsize=(6,6))
    w = [0.1, 1.0, 10.0]
    for nw in range(len(w)):
        ax.plot(evolve(W=w[nw])[-1], c=f'C{nw}', label=f'W={w[nw]:.1f}')
        ax.plot(X, analytic_sol(X,Mx//2,W=w[nw]), c=f'C{nw}', ls='--',label=f'W={w[nw]:.1f} solution')
    ax.set_xlabel('x',fontsize=12)
    ax.set_ylabel(r'$\eta$',fontsize=12)
    plt.legend()
    plt.savefig('W.png')
    plt.close()

    # kappa dependence
    fig, ax = plt.subplots(1, 1, figsize=(6,6))
    k = [0.1, 1.0, 10.0]
    for nk in range(len(k)):
        ax.plot(evolve(kappa=k[nk])[-1], c=f'C{nk}', label=f'kappa={k[nk]:.1f}')
        ax.plot(X, analytic_sol(X,Mx//2,kappa=k[nk]), c=f'C{nk}', ls='--',label=f'kappa={k[nk]:.1f} solution')
    ax.set_xlabel('x',fontsize=12)
    ax.set_ylabel(r'$\eta$',fontsize=12)
    plt.legend()
    plt.savefig('kappa.png')
