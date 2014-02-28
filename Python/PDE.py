from __future__ import division
from pylab import *

Nr = 1000; Nz = 10
Rmax = 2.5; Zmax = 1E-5
dR = Rmax / Nr; dZ = Zmax / Nz
R = np.linspace(0.1, Rmax, Nr)

a = -2j/dZ
b = 1/(2*dR)



def createA(A, E):
    for i in range(0, E.size):
        rI = abs(R[i])
        A[i, i] = rI *(a + 2*(b**2) - E[i]*np.conj(E[i]))
        x = b
        if i >= 1: A[i-1, i] = x
        if i < E.size-1: A[i, i+1] = -x
        if i >= 2: A[i-2, i] = -1*rI *(b**2)
        if i < E.size-2: A[i+2, i] = -1*rI *(b**2)


    A[0, 0] = A[0, 1] = 0.5*(a + (b**2) - E[0]*np.conj(E[0]))
    A[0, 2] = -(b**2)
    #
    A[1, 0] = 0
    A[1, 1] = a - (b**2) - E[1]*np.conj(E[1]) - b / R[1]
    A[1, 2] = 2*b**2 - b / R[1]
    A[1, 3] = -b**2
    A[1, 1] = a + b**2 - E[1]*np.conj(E[1]);

def createB(B, E):
    for i in range(0, E.size):
        rI = abs(R[i])
        B[i, i] = rI *(a - 2*(b**2) + E[i]*np.conj(E[i]))
        x = b
        if i >= 1: B[i-1, i] = -x
        if i < E.size-1: B[i, i+1] = x
        if i >= 2: B[i-2, i] = rI *b**2
        if i < E.size-2: B[i+2, i] = rI *b**2

    B[0, 0] = B[0, 1] = 0.5*(a - (b**2) + E[0]*np.conj(E[0]))
    B[0, 2] = (b**2)

    B[1, 0] = 0
    B[1, 1] = a + (b**2) + E[1]*np.conj(E[1]) + b / R[1]
    B[1, 2] = -2*b**2 - b / R[1]
    B[1, 3] = b**2
    B[1, 1] = a - b**2 + E[1]*np.conj(E[1])

def pde_solve():

    E=np.zeros((Nz+1, Nr),dtype='complex')
    G = np.exp(-np.square(R))
    E[0, 1:-1] = G[0:-2]
    E[0, 0] = G[0]

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


    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    # ax1.set_ylim([0, 1.1])
    # ax1.set_xlim([0, 2])

    plt.xlabel("r")
    plt.ylabel("Intensity")
# ax1.semilogx(data[:,1],data[:,2])
#     plt.arrow(0.05, 1.01, 0.55*R[-1], 0, ax1)

    ax2 = plt.axes([.65, .55, .2, .3])
# ax2.semilogx(data[3:8,1],data[3:8,2])
    plt.setp(ax2, xticks=[], yticks=[])
    ax1.plot(R, E0)
    ax1.plot(R, E1)
    ax1.plot(R, EN)

    Nm = 10
    low = 0#Nr / 2 - Nm
    high = Nm#Nr / 2 + Nm
    ax2.plot(R[low:high], E0[low:high])
    ax2.plot(R[low:high], E1[low:high])
    ax2.plot(R[low:high], EN[low:high])


    print(np.trapz(E0))
    print(np.trapz(EN))
    # plot(R, E[2, :])
    # plot(R, E[3, :])
    # plot(R, E[3, :])
    # plot(R, E[1, :])
    # plot(R, E[-1, :])
    show()

pde_solve()