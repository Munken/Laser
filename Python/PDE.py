from __future__ import division
from pylab import *

Nr = 2000; Nz = 50
Rmax = 2.5; Zmax = 0.0001
dR = Rmax / Nr; dZ = Zmax / Nz
R = np.linspace(1E-15, Rmax, Nr)

a = -2j/dZ
b = 1/(2*dR)

def createA(A, E):
    for i in range(0, E.size):
        A[i, i] = a + 2*(b**2) - E[i]*np.conj(E[i])
        x = b/R[i]
        if i >= 1: A[i-1, i] = -x
        if i < E.size-1: A[i, i+1] = x
        if i >= 2: A[i-2, i] = A[i, i-2] = -b**2


    A[0, 1] = 0
    A[0, 2] = -2*(b**2)
    A[1, 1] = a + b**2 - E[1]*np.conj(E[1]);

def createB(B, E):
    for i in range(0, E.size):
        B[i, i] = a - 2*(b**2) + E[i]*np.conj(E[i])
        x = b/R[i]
        if i >= 1: B[i-1, i] = x
        if i < E.size-1: B[i, i+1] = -x
        if i >= 2: B[i-2, i] = B[i, i-2] = b**2

    B[0, 1] = 0
    B[0, 2] = 2*(b**2)
    B[1, 1] = a - b**2 + E[1]*np.conj(E[1])

def pde_solve():


    E=np.zeros((Nz, Nr),dtype='complex')
    E[0, :] = np.exp(-np.square(R))

    A = np.zeros((Nr, Nr),dtype='complex')
    B = np.zeros((Nr, Nr),dtype='complex')
    for i in range(1, Nz):
        X = E[i-1, :]
        createA(A, X)
        createB(B, X)

        BE = np.dot(B, X)

        E[i, :] = np.linalg.solve(A, BE)


    # E[:] = E[:] * np.conj(E[:])
    E0 = E[0, :]
    E1 = E[1, :]
    EN = E[-1, :]

    E0 = E0 * np.conj(E0)
    E1 = E1 * np.conj(E1)
    EN = EN * np.conj(EN)
    plot(R, E0)
    plot(R, E1)
    plot(R, EN)
    # plot(R, E[2, :])
    # plot(R, E[3, :])
    # plot(R, E[3, :])
    # plot(R, E[1, :])
    # plot(R, E[-1, :])
    show()

pde_solve()