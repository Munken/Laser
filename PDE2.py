from pylab import *

Nr = 500; Nz = 1000
Rmax = 2.5; Zmax = 0.001
dR = Rmax / Nr; dZ = Zmax / Nz
R = np.linspace(1E-3, Rmax, Nr)
R2 = dR ** 2

def solve_PDE2():
    E=np.zeros((Nz+1, Nr),dtype='complex')
    E[0,:] = 1E-3*np.exp(-np.square(R))


    n = 1
    for i in range(1, Nr-1):
        E[n,i] = 1j*dZ*((E[n-1, i+1] - 2*E[n-1, i] + E[n-1, i-1]) / R2 + R[i]*(E[n-1, i+1] - E[n-1, i-1]/dR) + np.abs(E[n-1,i])**2*E[n-1,i]) + E[n-1, i]

    E[n,0] = 1j*dZ*((E[n-1,2]+E[n-1,1]-2*E[n-1,0])/ R2 + np.abs(E[n-1,0])**2*E[n-1,0]) + E[n-1,0]
    E[n,Nr-1] = 1j*dZ*((E[n-1,-3]+E[n-1,-1]-2*E[n-1,-2])/ R2 + R[-1]*(E[n-1, -2] - E[n-1, -1]/dR) + (np.abs(E[n-1,-1])**2)*E[n-1,-1]) + E[n-1,-1]

    for n in range(2, Nz+1):
        for i in range(1, Nr-1):
             E[n,i] = (2j*dZ*((E[n-1, i+1] - 2*E[n-1, i] + E[n-1, i-1]) / R2 + R[i]*(E[n-1, i+1] - E[n-1, i-1]/dR) + (np.abs(E[n-1,i])**2)*E[n-1,i]) + 4*E[n-1, i] - E[n-2, i])/3

        E[n,0] = (2j*dZ*((E[n-1,2]+E[n-1,1]-2*E[n-1,0])/ R2 + np.abs(E[n-1,0])**2*E[n-1,0]) + 4*E[n-1, 0] - E[n-2, 0])/3
        E[n,Nr-1] = (2j*dZ*((E[n-1,-3]+E[n-1,-1]-2*E[n-1,-2])/ R2 + R[-1]*(E[n-1, -2] - E[n-1, -1]/dR) + (np.abs(E[n-1,-1])**2)*E[n-1,-1]) + 4*E[n-1, -1] - E[n-2, -1])/3


    plot(R, np.abs(E[0, :]))
    # plot(R, E[1, :])
    # plot(R, E[200, :])
    # plot(R, E[300, :])
    # plot(R, E[400, :])
    plot(R, np.abs(E[-1, :]))




    print(np.trapz(np.abs(E[0, :])))
    print(np.trapz(np.abs(E[-1, :])))
    # plot(R, E[-1, :])
    show()



solve_PDE2()