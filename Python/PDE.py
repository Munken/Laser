from __future__ import division
from pylab import *

Nr = 500; Nz = 10
Rmax = 1; Zmax = 1E-5
dR = Rmax / Nr; dZ = Zmax / Nz
R = np.linspace(-Rmax, Rmax, Nr)

a = -2j/dZ
b = 1/(2*dR)



def createA(A, E):
    for i in range(0, E.size):
        rI = abs(R[i])
        A[i, i] = rI *(a + 2*(b**2) - E[i]*np.conj(E[i]))
        x = b
        if i >= 1: A[i-1, i] = x
        if i < E.size-1: A[i, i+1] = -x
        if i >= 2: A[i-2, i] = rI *(-b**2)
        if i < E.size-2: A[i+2, i] = rI *(-b**2)


    # A[0, 0] = A[0, 1] = 0.5*(a + (b**2) - E[0]*np.conj(E[0]))
    # A[0, 2] = -(b**2)
    #
    # A[1, 0] = 0
    # A[1, 1] = a - (b**2) - E[1]*np.conj(E[1]) - b / R[1]
    # A[1, 2] = 2*b**2 - b / R[1]
    # A[1, 3] = -b**2
    # A[1, 1] = a + b**2 - E[1]*np.conj(E[1]);

def createB(B, E):
    for i in range(0, E.size):
        rI = abs(R[i])
        B[i, i] = rI *(a - 2*(b**2) + E[i]*np.conj(E[i]))
        x = b
        if i >= 1: B[i-1, i] = -x
        if i < E.size-1: B[i, i+1] = x
        if i >= 2: B[i-2, i] = rI *b**2
        if i < E.size-2: B[i+2, i] = rI *b**2

    # B[0, 0] = B[0, 1] = 0.5*(a - (b**2) + E[0]*np.conj(E[0]))
    # B[0, 2] = (b**2)
    #
    # B[1, 0] = 0
    # B[1, 1] = a + (b**2) + E[1]*np.conj(E[1]) + b / R[1]
    # B[1, 2] = -2*b**2 - b / R[1]
    # B[1, 3] = b**2
    # B[1, 1] = a - b**2 + E[1]*np.conj(E[1])

def pde_solve():

    E=np.zeros((Nz+1, Nr),dtype='complex')
    E[0, :] = np.exp(-np.square(R)*20)

    A = np.zeros((Nr, Nr),dtype='complex')
    B = np.zeros((Nr, Nr),dtype='complex')
    for i in range(1, Nz+1):
        X = E[i-1, :]
        createA(A, X)
        createB(B, X)

        BE = np.dot(B, X)

        E[i, :] = np.linalg.solve(A, BE)

        if i % 200 == 0:
            Ei = E[i,:]
            Ei = Ei * np.conj(Ei)
            print("i ", i, " trapz", np.trapz(Ei))


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

    print(np.trapz(E0))
    print(np.trapz(EN))
    # plot(R, E[2, :])
    # plot(R, E[3, :])
    # plot(R, E[3, :])
    # plot(R, E[1, :])
    # plot(R, E[-1, :])
    show()

pde_solve()